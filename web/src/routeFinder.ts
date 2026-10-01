import { hohmann, route_xy, type Planet } from "./engine/engine.js";
import { COLOURS } from "./planets.js";
import type { Item } from "./plot/items.js";
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

const LEG_COLOURS = ["#b07cff", "#5b8def", "#2bb5a8", "#5fd39a", "#c6e05a", "#f0b85a"];

/** Route-finder panel: runs searches in a worker, lists results, and draws the chosen route. */
export class RouteFinder {
  /** Called with plot items and the date to show planets at, when a route is chosen or cleared. */
  onShow: (items: Item[], day: number | null) => void = () => {};

  private worker: Worker | null = null;
  private routes: FoundRoute[] = [];
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

    for (const s of [this.from, this.to]) s.innerHTML = bodies.map((b, k) => `<option value="${k}">${b.name}</option>`).join("");
    this.from.value = "2";
    this.to.value = "5";
    this.from.addEventListener("change", () => (this.defaults(), this.checkSequence()));
    this.to.addEventListener("change", () => (this.defaults(), this.checkSequence()));
    this.sequence.addEventListener("input", () => this.checkSequence());
    this.button.addEventListener("click", () => (this.worker ? this.stop("Stopped.") : this.search()));
    this.list.addEventListener("click", (e) => {
      const row = (e.target as HTMLElement).closest<HTMLElement>("[data-k]");
      if (row) this.show(Number(row.dataset.k));
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
    this.worker = new Worker(new URL("./routeWorker.ts", import.meta.url), { type: "module" });
    this.worker.onmessage = (e: MessageEvent<WorkerMessage>) => {
      const m = e.data;
      if (m.type === "progress") {
        this.status.textContent = m.status;
      } else {
        const what = flybys.length ? `launch windows for ${this.letters([from, ...flybys, to])}` : "routes";
        this.stop(`${m.routes.length} ${what} in ${(m.ms / 1000).toFixed(1)} s. ${this.status.textContent}`);
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

  private setRoutes(routes: FoundRoute[]) {
    this.routes = routes;
    this.list.innerHTML = routes.length
      ? `<table><thead><tr><th>Sequence</th><th>Launch</th><th class="num">Years</th><th class="num">Δv m/s</th></tr></thead><tbody>${routes
          .map((r, k) => {
            const years = (r.dates[r.dates.length - 1] - r.dates[0]) / 365.25;
            return `<tr data-k="${k}"><td>${this.letters(r.bodies)}</td><td>${fmtDate(r.dates[0])}</td><td class="num">${years.toFixed(1)}</td><td class="num">${r.dv.toFixed(0)}</td></tr>`;
          })
          .join("")}</tbody></table>`
      : "No route found. Try a wider launch window, a later arrival, or more flybys.";
    if (routes.length) this.show(0);
    else {
      this.events.textContent = "";
      this.onShow([], null);
    }
  }

  private show(k: number) {
    const r = this.routes[k];
    this.list.querySelectorAll<HTMLElement>("tr[data-k]").forEach((tr) => tr.classList.toggle("selected", tr.dataset.k === String(k)));
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
    this.events.textContent = `${lines.join("\n")}\n\nTotal Δv ${r.dv.toFixed(0)} m/s · ${years.toFixed(2)} years\nPlanets shown at launch.`;
  }
}
