import { hohmann, route_xy, type Planet } from "./engine/engine.js";
import { COLOURS } from "./planets.js";
import type { Item } from "./plot/items.js";
import { Refinement, type RefinedTrajectory } from "./refine.js";
import { parseFlybys } from "./sequence.js";
import { daysOf, fmtDate, parseDate } from "./time.js";

export interface SearchRequest {
  from: number;
  to: number;
  launchStart: number;
  launchEnd: number;
  maxFlybys: number;
  capture: boolean;
  /** m/s per year of flight. */
  timeCost: number;
  maxTof: number;
  step: number;
  margin: number;
  /** Planets to fly past, or empty to search for the best sequence. */
  flybys: number[];
}

/** A route as the worker returns it. Costs in m/s, dates in days after J2000. */
export interface FoundRoute {
  bodies: number[];
  dates: number[];
  c3: number;
  departure: number;
  flybys: number[];
  arrival: number;
  dv: number;
  score: number;
}

export type WorkerMessage =
  | { type: "progress"; status: string }
  | { type: "done"; routes: FoundRoute[]; ms: number };

/** A listed route: the grid search's version and, once SADE has run on it, the refined one. */
interface Entry {
  id: number;
  route: FoundRoute;
  refined: RefinedTrajectory | null;
  state: "queued" | "running" | "done" | "stopped" | null;
}

const LEG_COLOURS = ["#b07cff", "#5b8def", "#2bb5a8", "#5fd39a", "#c6e05a", "#f0b85a"];

/** Δv to rank by: the refined trajectory's once there is one, else the grid route's. */
const dvOf = (e: Entry) => e.refined?.dv ?? e.route.dv;

/**
 * Route-finder panel: runs searches in a worker, lists the routes, then refines each with SADE
 * in turn (best first), re-sorting the list by Δv as refinements finish. Draws the chosen route.
 */
export class RouteFinder {
  /** Called with plot items and the date to show planets at, when a route is chosen or cleared. */
  onShow: (items: Item[], day: number | null) => void = () => {};

  private worker: Worker | null = null;
  /** The search the listed routes came from. */
  private searched: SearchRequest | null = null;
  /** Listed routes, in display order. */
  private entries: Entry[] = [];
  private selectedId = -1;
  /** Set when the user picks a row other than the one being refined: the plot then stays on it. */
  private pinned = false;
  private queue: Entry[] = [];
  private running: { entry: Entry; refinement: Refinement } | null = null;
  private queueSize = 0;

  private readonly from: HTMLSelectElement;
  private readonly to: HTMLSelectElement;
  private readonly launchStart: HTMLInputElement;
  private readonly launchEnd: HTMLInputElement;
  private readonly arriveBy: HTMLInputElement;
  private readonly maxFlybys: HTMLInputElement;
  private readonly sequence: HTMLInputElement;
  private readonly sequenceHint: HTMLElement;
  private readonly arrival: HTMLSelectElement;
  private readonly timeCost: HTMLInputElement;
  private readonly step: HTMLInputElement;
  private readonly button: HTMLButtonElement;
  private readonly status: HTMLElement;
  private readonly list: HTMLElement;
  private readonly events: HTMLElement;
  private readonly refinePanel: HTMLElement;
  private readonly islands: HTMLInputElement;
  private readonly rounds: HTMLInputElement;
  private readonly refineButton: HTMLButtonElement;
  private readonly refineStatus: HTMLElement;
  private readonly copyButton: HTMLButtonElement;

  constructor(root: HTMLElement, private readonly bodies: Planet[]) {
    const $ = <T extends HTMLElement>(id: string) => root.querySelector<T>(`#${id}`)!;
    this.from = $("rf-from");
    this.to = $("rf-to");
    this.launchStart = $("rf-launch-start");
    this.launchEnd = $("rf-launch-end");
    this.arriveBy = $("rf-arrive-by");
    this.maxFlybys = $("rf-flybys");
    this.sequence = $("rf-sequence");
    this.sequenceHint = $("rf-sequence-hint");
    this.arrival = $("rf-arrival");
    this.timeCost = $("rf-time-cost");
    this.step = $("rf-step");
    this.button = $("rf-search");
    this.status = $("rf-status");
    this.list = $("rf-results");
    this.events = $("rf-events");
    this.refinePanel = $("rf-refine");
    this.islands = $("rf-islands");
    this.rounds = $("rf-rounds");
    this.refineButton = $("rf-refine-go");
    this.refineStatus = $("rf-refine-status");
    this.copyButton = $("rf-copy");

    for (const s of [this.from, this.to]) s.innerHTML = bodies.map((b, k) => `<option value="${k}">${b.name}</option>`).join("");
    this.from.value = "2";
    this.to.value = "5";
    this.islands.value = String(Math.max(1, Math.min(16, (navigator.hardwareConcurrency || 4) - 1)));
    this.from.addEventListener("change", () => (this.defaults(), this.checkSequence()));
    this.to.addEventListener("change", () => (this.defaults(), this.checkSequence()));
    this.sequence.addEventListener("input", () => this.checkSequence());
    this.button.addEventListener("click", () => (this.worker ? this.stop("Stopped.") : this.search()));
    this.refineButton.addEventListener("click", () => (this.running ? this.stopQueue("Stopped.") : this.refineAll()));
    this.copyButton.addEventListener("click", () => {
      const e = this.entry(this.selectedId);
      if (e?.refined) void navigator.clipboard.writeText(JSON.stringify(e.refined.chromosome));
    });
    this.list.addEventListener("click", (e) => {
      const row = (e.target as HTMLElement).closest<HTMLElement>("[data-id]");
      if (!row) return;
      const id = Number(row.dataset.id);
      this.pinned = !!this.running && id !== this.running.entry.id;
      this.select(id);
    });
    this.defaults();
  }

  /** Parses the sequence field, flagging errors and showing what it was read as. */
  private checkSequence(): number[] | null {
    const names = this.bodies.map((b) => b.name);
    const from = Number(this.from.value);
    const to = Number(this.to.value);
    const parsed = parseFlybys(this.sequence.value, names, from, to);
    const bad = typeof parsed === "string";
    this.sequence.setAttribute("aria-invalid", String(bad));
    this.maxFlybys.disabled = !bad && parsed.length > 0;
    this.sequenceHint.textContent = bad
      ? parsed
      : parsed.length
        ? `Route: ${[from, ...parsed, to].map((k) => names[k]).join(" → ")}. Gives the best route in each launch window.`
        : "Optional: the planets to fly past, e.g. VVEJ or EVVEJS (M is Mars; write Me for Mercury). Gives the best route in each launch window.";
    return bad ? null : parsed;
  }

  /** Sensible ranges for the chosen pair: a launch window of 2–5 years, arrivals up to about twice
   *  the Hohmann time later, and a grid step that keeps the search quick. */
  private defaults() {
    const a = this.bodies[Number(this.from.value)];
    const b = this.bodies[Number(this.to.value)];
    const t = hohmann(a.a_au, b.a_au);
    const th = t.tof_days;
    t.free();
    const synodic = a === b ? 365 : 1 / Math.abs(1 / a.period_days - 1 / b.period_days);
    const start = daysOf(new Date());
    const window = Math.min(Math.max(synodic, 730), 1826);
    this.launchStart.value = fmtDate(start);
    this.launchEnd.value = fmtDate(start + window);
    this.arriveBy.value = fmtDate(start + window + Math.min(2 * th + 365, 15 * 365.25));
    this.step.value = String(Math.min(Math.max(Math.round(th / 150), 3), 15));
  }

  private search() {
    const from = Number(this.from.value);
    const to = Number(this.to.value);
    const launchStart = parseDate(this.launchStart.value);
    const launchEnd = parseDate(this.launchEnd.value);
    const arriveBy = parseDate(this.arriveBy.value);
    const flybys = this.checkSequence();
    if (!flybys) {
      this.status.textContent = "Fix the sequence first.";
      return;
    }
    const req: SearchRequest = {
      from, to, launchStart, launchEnd,
      maxFlybys: this.maxFlybys.valueAsNumber,
      capture: this.arrival.value === "capture",
      timeCost: this.timeCost.valueAsNumber || 0,
      maxTof: arriveBy - launchEnd,
      step: this.step.valueAsNumber,
      margin: 0.15,
      flybys,
    };
    if (from === to || !(launchEnd > launchStart && req.maxTof > 0 && req.step > 0 && req.maxFlybys >= 0)) {
      this.status.textContent = "Check the planets, dates and settings.";
      return;
    }
    this.stopQueue("");
    this.worker = new Worker(new URL("./routeWorker.ts", import.meta.url), { type: "module" });
    this.worker.onmessage = (e: MessageEvent<WorkerMessage>) => {
      const m = e.data;
      if (m.type === "progress") {
        this.status.textContent = m.status;
      } else {
        const what = flybys.length ? `launch windows for ${this.letters([from, ...flybys, to])}` : "routes";
        this.stop(`${m.routes.length} ${what} in ${(m.ms / 1000).toFixed(1)} s. ${this.status.textContent}`);
        this.searched = req;
        this.setRoutes(m.routes);
      }
    };
    this.worker.postMessage(req);
    this.button.textContent = "Stop";
    this.status.textContent = "Searching…";
  }

  private stop(status: string) {
    this.worker?.terminate();
    this.worker = null;
    this.button.textContent = "Find routes";
    this.status.textContent = status;
  }

  /** Compact sequence, e.g. EVVEJS (Mercury is Me, as in the sequence field). */
  private letters(bodies: number[]) {
    return bodies.map((b) => (this.bodies[b].name === "Mercury" ? "Me" : this.bodies[b].name[0])).join("");
  }

  private entry(id: number) {
    return this.entries.find((e) => e.id === id);
  }

  /** Lists new search results and starts refining them, best first. */
  private setRoutes(routes: FoundRoute[]) {
    this.entries = routes.map((route, id) => ({ id, route, refined: null, state: null }));
    this.pinned = false;
    this.refinePanel.hidden = !routes.length;
    this.refineStatus.textContent = "";
    if (!routes.length) {
      this.list.textContent = "No route found. Try a wider launch window, a later arrival, or more flybys.";
      this.events.textContent = "";
      this.copyButton.hidden = true;
      this.onShow([], null);
      return;
    }
    this.select(0);
    this.refineAll();
  }

  /** Text of a row's changing cells: launch, years, and SADE Δv. */
  private cells(e: Entry): [string, string, string] {
    const t = e.refined ?? e.route;
    const years = (t.dates[t.dates.length - 1] - t.dates[0]) / 365.25;
    const sade =
      e.state === "queued" ? "queued"
      : !e.refined ? "–"
      : `${e.refined.dv.toFixed(0)}${e.state === "running" ? " …" : e.state === "stopped" ? " *" : ""}`;
    return [fmtDate(t.dates[0]), years.toFixed(1), sade];
  }

  /** Updates one row in place, so rows being clicked aren't replaced under the pointer. */
  private updateRow(e: Entry) {
    const tr = this.list.querySelector<HTMLTableRowElement>(`tr[data-id="${e.id}"]`);
    if (!tr) return;
    const [launch, years, sade] = this.cells(e);
    tr.cells[1].textContent = launch;
    tr.cells[2].textContent = years;
    tr.cells[4].textContent = sade;
  }

  private renderTable() {
    const rows = this.entries.map((e) => {
      const [launch, years, sade] = this.cells(e);
      const cls = [e.id === this.selectedId ? "selected" : "", e.state === "running" ? "running" : ""].join(" ").trim();
      const sadeCls = e.state === "running" ? "num live" : e.refined ? "num" : "num dim";
      return (
        `<tr data-id="${e.id}" class="${cls}"><td>${this.letters(e.route.bodies)}</td><td>${launch}</td>` +
        `<td class="num">${years}</td><td class="num">${e.route.dv.toFixed(0)}</td><td class="${sadeCls}">${sade}</td></tr>`
      );
    });
    this.list.innerHTML =
      `<table><thead><tr><th>Sequence</th><th>Launch</th><th class="num">Years</th>` +
      `<th class="num" title="Flyby-only grid route, powered flybys">Grid Δv</th>` +
      `<th class="num" title="After SADE on the deep-space-burn model">SADE Δv</th></tr></thead>` +
      `<tbody>${rows.join("")}</tbody></table>`;
  }

  private select(id: number) {
    this.selectedId = id;
    this.renderTable();
    const e = this.entry(id);
    if (!e) return;
    if (e.refined) this.showRefined(e.refined, e.route);
    else this.showGrid(e.route);
  }

  /** Queues every route not yet refined, in list order, and starts the first. */
  private refineAll() {
    this.queue = this.entries.filter((e) => e.state !== "done");
    for (const e of this.queue) e.state = "queued";
    this.queueSize = this.queue.length;
    this.refineButton.textContent = "Stop";
    this.renderTable();
    this.next();
  }

  /** Refines the next queued route with seeded SADE on the deep-space-burn model, one island per
   *  worker, updating its row and the plot whenever any island improves on it. */
  private next() {
    const q = this.searched;
    const entry = this.queue.shift();
    if (!q || !entry) {
      this.running = null;
      this.refineButton.textContent = "Refine all";
      if (q && this.queueSize) this.refineStatus.textContent = `Refined ${this.queueSize} routes; list sorted by Δv.`;
      // Finish on the new best, unless the user has picked a route to look at.
      if (!this.pinned && this.entries.length) this.select(this.entries[0].id);
      return;
    }
    const r = entry.route;
    const rounds = Math.max(1, this.rounds.valueAsNumber || 1);
    const islands = Math.max(1, Math.min(32, this.islands.valueAsNumber || 1));
    const refinement = new Refinement(
      {
        rounds, bodies: r.bodies, dates: r.dates, launchStart: q.launchStart, launchEnd: q.launchEnd,
        capture: q.capture, timeCost: q.timeCost, maxTof: q.maxTof, step: q.step,
      },
      islands,
    );
    this.running = { entry, refinement };
    entry.state = "running";
    if (!this.pinned) this.selectedId = entry.id;
    this.renderTable();

    const label = `${this.queueSize - this.queue.length} of ${this.queueSize}: ${this.letters(r.bodies)} ${fmtDate(r.dates[0])}`;
    refinement.onBest = (t) => {
      entry.refined = t;
      this.updateRow(entry);
      if (this.selectedId === entry.id) this.showRefined(t, r);
    };
    refinement.onProgress = (p) => {
      const n = (v: number) => Math.round(v).toLocaleString();
      this.refineStatus.textContent =
        `Refining ${label} · round ${Math.min(p.rounds + 1, rounds)} of ${rounds} · ${p.islands} islands · ` +
        `${(p.evaluations / 1e6).toFixed(2)} M evaluations · ${p.seconds.toFixed(0)} s\n` +
        `Score: grid ${n(r.score)} → seed ${n(p.seedScore)} → best ${n(p.bestScore)} m/s`;
    };
    refinement.onDone = () => {
      entry.state = "done";
      this.sortByDv();
      this.next();
    };
  }

  /** Stops refining; the route in progress keeps its best result so far, marked as stopped. */
  private stopQueue(status: string) {
    if (this.running) {
      this.running.refinement.stop();
      this.running.entry.state = this.running.entry.refined ? "stopped" : null;
      this.running = null;
    }
    for (const e of this.queue) e.state = null;
    this.queue = [];
    this.refineButton.textContent = "Refine all";
    this.refineStatus.textContent = status ? `${status} ${this.refineStatus.textContent}` : "";
    this.renderTable();
  }

  /** Stable re-sort by Δv, refined where available. */
  private sortByDv() {
    this.entries = this.entries
      .map((e, k) => ({ e, k }))
      .sort((a, b) => dvOf(a.e) - dvOf(b.e) || a.k - b.k)
      .map(({ e }) => e);
    this.renderTable();
  }

  /** Draws a grid (flyby-only) route and lists its events. */
  private showGrid(r: FoundRoute) {
    this.copyButton.hidden = true;
    const n = 120;
    const xy = route_xy(new Uint32Array(r.bodies), new Float64Array(r.dates), n);
    const legs = r.bodies.length - 1;
    const per = 2 * (n + 1);
    const items: Item[] = [];
    for (let l = 0; l < legs; l++) {
      items.push({ kind: "line", xy: xy.subarray(per * l, per * (l + 1)), color: LEG_COLOURS[l % LEG_COLOURS.length], width: 2 });
    }
    r.bodies.forEach((b, e) => {
      const p = e < legs ? per * e : per * legs - 2;
      const shape = e === 0 ? "triangle" : e === legs ? "square" : "circle";
      items.push({ kind: "marker", x: xy[p], y: xy[p + 1], shape, color: COLOURS[this.bodies[b].name], size: 6, tag: String(e + 1) });
    });
    this.onShow(items, r.dates[0]);

    const name = (b: number) => this.bodies[b].name.padEnd(8);
    const lines = r.bodies.map((b, e) => {
      const date = fmtDate(r.dates[e]);
      if (e === 0) return `${e + 1}  ${date}  launch  ${name(b)} C3 ${r.c3.toFixed(1)} km²/s², ${r.departure.toFixed(0)} m/s`;
      if (e === legs) return `${e + 1}  ${date}  ${r.arrival > 0 ? "capture" : "arrive "} ${name(b)} ${r.arrival.toFixed(0)} m/s`;
      return `${e + 1}  ${date}  flyby   ${name(b)} ${r.flybys[e - 1].toFixed(0)} m/s`;
    });
    const years = (r.dates[legs] - r.dates[0]) / 365.25;
    this.events.textContent = `Grid route (powered flybys, not yet refined)\n${lines.join("\n")}\n\nTotal Δv ${r.dv.toFixed(0)} m/s · ${years.toFixed(2)} years\nPlanets shown at launch.`;
  }

  /** Draws a deep-space-burn trajectory and lists its events in time order, numbered as on the plot. */
  private showRefined(t: RefinedTrajectory, grid: FoundRoute) {
    this.copyButton.hidden = false;
    const legs = t.bodies.length - 1;
    const per = t.xy.length / legs;
    const items: Item[] = [];
    for (let l = 0; l < legs; l++) {
      items.push({ kind: "line", xy: t.xy.subarray(per * l, per * (l + 1)), color: LEG_COLOURS[l % LEG_COLOURS.length], width: 2 });
    }
    const name = (b: number) => this.bodies[b].name;
    const lines: string[] = [];
    let tag = 1;
    const mark = (x: number, y: number, shape: "triangle" | "square" | "circle" | "star" | "cross", color: string) =>
      items.push({ kind: "marker", x, y, shape, color, size: shape === "cross" ? 4 : 6, tag: String(tag) });
    const line = (date: number, what: string, body: string, detail: string) =>
      lines.push(`${String(tag++).padStart(2)}  ${fmtDate(date)}  ${what.padEnd(7)} ${body.padEnd(17)} ${detail}`);

    mark(t.xy[0], t.xy[1], "triangle", COLOURS[name(t.bodies[0])]);
    line(t.dates[0], "launch", name(t.bodies[0]), `C3 ${t.c3.toFixed(1)} km²/s², ${t.departure.toFixed(0)} m/s`);
    for (let l = 0; l < legs; l++) {
      const big = t.dsms[l] >= 1;
      mark(t.dsmXy[2 * l], t.dsmXy[2 * l + 1], big ? "star" : "cross", big ? "#ff4d4d" : "#8a94a6");
      line(t.dsmDates[l], "burn", `${name(t.bodies[l])}→${name(t.bodies[l + 1])}`, big ? `${t.dsms[l].toFixed(0)} m/s` : "negligible");
      const end = per * (l + 1) - 2;
      const b = t.bodies[l + 1];
      if (l + 1 < legs) {
        mark(t.xy[end], t.xy[end + 1], "circle", COLOURS[name(b)]);
        line(t.dates[l + 1], "flyby", name(b), "unpowered");
      } else {
        mark(t.xy[end], t.xy[end + 1], "square", COLOURS[name(b)]);
        line(t.dates[l + 1], t.arrival > 0 ? "capture" : "arrive", name(b), t.arrival > 0 ? `${t.arrival.toFixed(0)} m/s` : "");
      }
    }
    this.onShow(items, t.dates[0]);

    const years = (t.dates[legs] - t.dates[0]) / 365.25;
    const burns = t.dsms.reduce((a, b) => a + b, 0);
    this.events.textContent =
      `Refined with SADE (deep-space burns, unpowered flybys)\n${lines.join("\n")}\n\n` +
      `Total Δv ${t.dv.toFixed(0)} m/s (deep-space burns ${burns.toFixed(0)}) · ${years.toFixed(2)} years\n` +
      `Grid route ${grid.dv.toFixed(0)} m/s with powered flybys. Star = burn, grey × = negligible (< 1 m/s).\n` +
      `Planets shown at launch.`;
  }
}
