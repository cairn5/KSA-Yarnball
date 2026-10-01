/** One SADE island: evolves a seeded population and streams its best trajectory as it improves. */
import init, { Refiner, type Trajectory } from "./engine/engine.js";
import type { IslandMessage, RefineRequest, RefinedTrajectory } from "./refine.js";

/** Points per coast and per arc when drawing. */
const DRAW_POINTS = 60;
/** Generations between progress reports. */
const CHUNK = 20;

const ready = init();
const send = (msg: IslandMessage, transfer: Transferable[] = []) => postMessage(msg, { transfer });

function plain(t: Trajectory): RefinedTrajectory {
  const out: RefinedTrajectory = {
    bodies: Array.from(t.bodies), dates: Array.from(t.dates), dsmDates: Array.from(t.dsm_dates),
    dsms: Array.from(t.dsms), c3: t.c3, departure: t.departure, arrival: t.arrival, dv: t.dv,
    score: t.score, xy: t.xy, dsmXy: t.dsm_xy, chromosome: Array.from(t.chromosome),
  };
  t.free();
  return out;
}

onmessage = async (e: MessageEvent<RefineRequest>) => {
  await ready;
  const q = e.data;
  const r = new Refiner(
    new Uint32Array(q.bodies), new Float64Array(q.dates), q.launchStart, q.launchEnd,
    q.capture, q.timeCost, q.maxTof, q.step, BigInt(q.seed),
  );
  const seed = r.describe(r.seed, DRAW_POINTS);
  if (seed) {
    const t = plain(seed);
    send({ type: "seed", island: q.island, trajectory: t }, [t.xy.buffer, t.dsmXy.buffer]);
  }
  let best = Infinity;
  while (r.rounds < q.rounds) {
    r.evolve(CHUNK);
    const score = r.best_score;
    if (score < best) {
      best = score;
      const d = r.describe(r.best, DRAW_POINTS);
      if (d) {
        const t = plain(d);
        send({ type: "best", island: q.island, trajectory: t }, [t.xy.buffer, t.dsmXy.buffer]);
      }
    }
    send({ type: "progress", island: q.island, generation: r.generation, evaluations: r.evaluations, rounds: r.rounds, best: score });
  }
  r.free();
  send({ type: "done", island: q.island });
};
