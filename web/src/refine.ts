/**
 * Short seeded SADE across islands: one Web Worker per island, each an independent population
 * started from the same route (pygmo's default archipelago is unconnected, so islands never
 * exchange solutions). The best trajectory across islands is reported whenever it improves.
 */

/** What an island needs: a route-finder route and the settings it was found with. */
export interface RefineRequest {
  island: number;
  seed: number;
  rounds: number;
  bodies: number[];
  dates: number[];
  launchStart: number;
  launchEnd: number;
  capture: boolean;
  timeCost: number;
  maxTof: number;
  step: number;
}

/** An MGA-1DSM trajectory. Costs in m/s, dates in days after J2000, drawing in AU. */
export interface RefinedTrajectory {
  bodies: number[];
  dates: number[];
  /** One burn per leg. */
  dsmDates: number[];
  dsms: number[];
  c3: number;
  departure: number;
  arrival: number;
  dv: number;
  score: number;
  /** Each leg's coast then arc, interleaved x, y: `xy.length / legs` values per leg. */
  xy: Float64Array;
  dsmXy: Float64Array;
  /** pykep mga_1dsm chromosome (direct encoding). */
  chromosome: number[];
}

export type IslandMessage =
  | { type: "seed"; island: number; trajectory: RefinedTrajectory }
  | { type: "best"; island: number; trajectory: RefinedTrajectory }
  | { type: "progress"; island: number; generation: number; evaluations: number; rounds: number; best: number }
  | { type: "done"; island: number };

export interface RefineProgress {
  islands: number;
  finished: number;
  /** Rounds completed by the slowest island. */
  rounds: number;
  evaluations: number;
  seedScore: number;
  bestScore: number;
  seconds: number;
}

export class Refinement {
  onBest: (t: RefinedTrajectory) => void = () => {};
  onProgress: (p: RefineProgress) => void = () => {};
  onDone: () => void = () => {};

  private readonly workers: Worker[] = [];
  private readonly stats: { evaluations: number; rounds: number; done: boolean }[] = [];
  private best: RefinedTrajectory | null = null;
  private stopped = false;
  private seedScore = Infinity;
  private readonly started = performance.now();

  constructor(base: Omit<RefineRequest, "island" | "seed">, islands: number) {
    for (let k = 0; k < islands; k++) {
      const w = new Worker(new URL("./refineWorker.ts", import.meta.url), { type: "module" });
      this.stats.push({ evaluations: 0, rounds: 0, done: false });
      w.onmessage = (e: MessageEvent<IslandMessage>) => this.receive(e.data);
      w.postMessage({ ...base, island: k, seed: 1 + k * 7919 } satisfies RefineRequest);
      this.workers.push(w);
    }
  }

  stop() {
    this.stopped = true;
    for (const w of this.workers) w.terminate();
  }

  private receive(m: IslandMessage) {
    // Messages already queued when the workers were stopped are dropped.
    if (this.stopped) return;
    const s = this.stats[m.island];
    if (m.type === "seed") {
      this.seedScore = m.trajectory.score;
      this.consider(m.trajectory);
    } else if (m.type === "best") {
      this.consider(m.trajectory);
    } else if (m.type === "progress") {
      s.evaluations = m.evaluations;
      s.rounds = m.rounds;
    } else {
      s.done = true;
    }
    this.onProgress({
      islands: this.stats.length,
      finished: this.stats.filter((x) => x.done).length,
      rounds: Math.min(...this.stats.map((x) => x.rounds)),
      evaluations: this.stats.reduce((a, x) => a + x.evaluations, 0),
      seedScore: this.seedScore,
      bestScore: this.best?.score ?? Infinity,
      seconds: (performance.now() - this.started) / 1000,
    });
    if (this.stats.every((x) => x.done)) {
      this.stop();
      this.onDone();
    }
  }

  private consider(t: RefinedTrajectory) {
    if (this.best && t.score >= this.best.score) return;
    this.best = t;
    this.onBest(t);
  }
}
