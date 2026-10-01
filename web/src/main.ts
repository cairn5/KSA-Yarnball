import init, { hohmann, lambert_leg, planets, porkchop } from "./engine/engine.js";
import { Plot } from "./plot/plot.js";
import { COLOURS, orbitLayer, planetMarkers, showElements } from "./planets.js";
import { PorkchopView, type Metric } from "./porkchop.js";
import type { Item } from "./plot/items.js";
import { RouteFinder } from "./routeFinder.js";
import { daysOf, fmtDate, parseDate } from "./time.js";

await init();

const $ = <T extends HTMLElement>(sel: string) => document.querySelector<T>(sel)!;

const plot = new Plot($<HTMLCanvasElement>("#plot"));
const bodies = planets();
const byName = (name: string) => bodies.find((b) => b.name === name)!;
plot.setLayer("orbits", orbitLayer(bodies));
plot.setLayer("transfer", []);
plot.setLayer("planets", planetMarkers(bodies, daysOf(new Date())));

const info = $<HTMLElement>("#info");
plot.onSelect = (id) => {
  const p = bodies.find((b) => b.name === id);
  info.hidden = !p;
  if (p) showElements(info, p);
};

// Porkchop panel.
const from = $<HTMLSelectElement>("#from");
const to = $<HTMLSelectElement>("#to");
const depStart = $<HTMLInputElement>("#dep-start");
const depDays = $<HTMLInputElement>("#dep-days");
const tofMin = $<HTMLInputElement>("#tof-min");
const tofMax = $<HTMLInputElement>("#tof-max");
const metric = $<HTMLSelectElement>("#metric");
const status = $<HTMLElement>("#status");
const legInfo = $<HTMLElement>("#leg");
const view = new PorkchopView($<HTMLCanvasElement>("#porkchop"), $<HTMLElement>("#readout"));

for (const select of [from, to]) select.innerHTML = bodies.map((b) => `<option>${b.name}</option>`).join("");
from.value = "Earth";
to.value = "Mars";
depStart.value = fmtDate(daysOf(new Date()));

/** Sensible ranges for a pair: one synodic period of departures, flight times around the Hohmann time. */
function defaultRanges() {
  const a = byName(from.value);
  const b = byName(to.value);
  const synodic = 1 / Math.abs(1 / a.period_days - 1 / b.period_days);
  const t = hohmann(a.a_au, b.a_au);
  depDays.value = String(Math.round(Math.min(Math.max(synodic, 200), 1000)));
  tofMin.value = String(Math.max(5, Math.round(0.4 * t.tof_days)));
  tofMax.value = String(Math.round(1.6 * t.tof_days));
  t.free();
}

function compute() {
  if (from.value === to.value) {
    status.textContent = "Pick two different planets.";
    return;
  }
  const dep0 = parseDate(depStart.value);
  const span = depDays.valueAsNumber;
  const t0 = tofMin.valueAsNumber;
  const t1 = tofMax.valueAsNumber;
  if (!(Number.isFinite(dep0) && span > 0 && t0 > 0 && t1 > t0)) {
    status.textContent = "Check the date and ranges.";
    return;
  }
  const [nDep, nTof] = [160, 120];
  const started = performance.now();
  const result = porkchop(byName(from.value), byName(to.value), dep0, dep0 + span, nDep, t0, t1, nTof);
  const grid = {
    depStart: dep0, depEnd: dep0 + span, nDep,
    tofMin: t0, tofMax: t1, nTof,
    c3: result.c3(), vinfArr: result.vinf_arr(),
  };
  result.free();
  status.textContent = `${(nDep * nTof).toLocaleString()} Lambert solves in ${(performance.now() - started).toFixed(0)} ms`;
  view.setGrid(grid);
  view.pickBest();
}

/** What each tab last drew: plot items and the date to show planets at. */
const shown: Record<string, [Item[], number | null]> = { porkchop: [[], null], routes: [[], null] };
let activeTab = "porkchop";

function show(tab: string, items: Item[], day: number | null) {
  shown[tab] = [items, day];
  if (tab !== activeTab) return;
  plot.setLayer("transfer", items);
  plot.setLayer("planets", planetMarkers(bodies, day ?? daysOf(new Date())));
}

view.onPick = (dep, tof) => {
  const a = byName(from.value);
  const b = byName(to.value);
  const leg = lambert_leg(a, b, dep, tof, 240);
  if (!leg) {
    show("porkchop", [], dep);
    legInfo.textContent = "No transfer for that cell.";
    return;
  }
  const { xy } = leg;
  const n = xy.length;
  show("porkchop", [
    { kind: "line", xy, color: "#ff5d8f", width: 2 },
    { kind: "marker", x: xy[0], y: xy[1], shape: "triangle", color: COLOURS[a.name], size: 6, tag: "1" },
    { kind: "marker", x: xy[n - 2], y: xy[n - 1], shape: "square", color: COLOURS[b.name], size: 5, tag: "2" },
  ], dep);
  legInfo.textContent =
    `1  ${fmtDate(dep)}  launch  ${a.name}   C3 ${leg.c3.toFixed(2)} km²/s²\n` +
    `2  ${fmtDate(dep + tof)}  arrive  ${b.name}   v∞ ${leg.vinf_arr.toFixed(2)} km/s\n` +
    `Planets shown at launch.`;
  leg.free();
};

from.addEventListener("change", () => (defaultRanges(), compute()));
to.addEventListener("change", () => (defaultRanges(), compute()));
for (const input of [depStart, depDays, tofMin, tofMax]) input.addEventListener("change", compute);
metric.addEventListener("change", () => {
  view.setMetric(metric.value as Metric);
  view.pickBest();
});

defaultRanges();
compute();

// Route finder panel, and switching between the two.
const finder = new RouteFinder($<HTMLElement>("#tab-routes"), bodies);
finder.onShow = (items, day) => show("routes", items, day);

for (const tab of document.querySelectorAll<HTMLButtonElement>("[data-tab]")) {
  tab.addEventListener("click", () => {
    activeTab = tab.dataset.tab!;
    for (const t of document.querySelectorAll<HTMLButtonElement>("[data-tab]")) {
      t.setAttribute("aria-selected", String(t === tab));
      $<HTMLElement>(`#tab-${t.dataset.tab}`).hidden = t !== tab;
    }
    show(activeTab, ...shown[activeTab]);
  });
}
