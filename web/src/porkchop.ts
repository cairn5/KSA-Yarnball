import { contours, niceLevels } from "./contours.js";
import { BG, FONT, INK } from "./plot/items.js";
import { dateOf, fmtDate } from "./time.js";

export type Metric = "c3" | "vinfArr" | "total";

/** A porkchop grid. Arrays are row-major: index = tofIndex * nDep + depIndex. Dates are days after J2000. */
export interface Grid {
  depStart: number;
  depEnd: number;
  nDep: number;
  tofMin: number;
  tofMax: number;
  nTof: number;
  c3: Float64Array;
  vinfArr: Float64Array;
}

const METRICS: Record<Metric, { label: string; unit: string; of: (c3: number, va: number) => number }> = {
  c3: { label: "Launch C3", unit: "km²/s²", of: (c3) => c3 },
  vinfArr: { label: "Arrival v∞", unit: "km/s", of: (_, va) => va },
  total: { label: "Total v∞", unit: "km/s", of: (c3, va) => Math.sqrt(c3) + va },
};

// Viridis, reversed: cheapest is brightest.
const STOPS = [
  [253, 231, 37],
  [94, 201, 98],
  [33, 145, 140],
  [59, 82, 139],
  [68, 1, 84],
];

function colour(t: number): [number, number, number] {
  const x = Math.min(Math.max(t, 0), 1) * (STOPS.length - 1);
  const i = Math.min(Math.floor(x), STOPS.length - 2);
  const f = x - i;
  const [a, b] = [STOPS[i], STOPS[i + 1]];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
}

const M = { left: 42, right: 8, top: 8, bottom: 30 };

/** Heatmap of departure date (x) against time of flight (y). Hover to read values, click to pick a transfer. */
export class PorkchopView {
  metric: Metric = "c3";
  /** Called with departure date and time of flight (days) when a cell is clicked. */
  onPick: (dep: number, tof: number) => void = () => {};

  private grid: Grid | null = null;
  private values = new Float64Array(0);
  private best = -1;
  private lines: { level: number; segs: number[] }[] = [];
  private picked: [number, number] | null = null;
  private readonly ctx: CanvasRenderingContext2D;
  private readonly image = document.createElement("canvas");

  constructor(private readonly canvas: HTMLCanvasElement, private readonly readout: HTMLElement) {
    this.ctx = canvas.getContext("2d")!;
    new ResizeObserver(() => this.draw()).observe(canvas);
    canvas.addEventListener("pointermove", (e) => this.hover(e.offsetX, e.offsetY));
    canvas.addEventListener("pointerleave", () => this.showCell(this.picked));
    canvas.addEventListener("click", (e) => {
      const cell = this.cellAt(e.offsetX, e.offsetY);
      if (cell) this.pick(...cell);
    });
  }

  setGrid(grid: Grid) {
    this.grid = grid;
    this.update();
  }

  setMetric(metric: Metric) {
    this.metric = metric;
    this.update();
  }

  /** Selects the cell with the lowest value of the current metric. */
  pickBest() {
    const g = this.grid;
    if (!g || this.best < 0) return;
    this.pick(this.depAt(this.best % g.nDep), this.tofAt(Math.floor(this.best / g.nDep)));
  }

  private pick(dep: number, tof: number) {
    this.picked = [dep, tof];
    this.draw();
    this.showCell(this.picked);
    this.onPick(dep, tof);
  }

  private depAt(i: number) {
    const g = this.grid!;
    return g.depStart + ((g.depEnd - g.depStart) * i) / Math.max(g.nDep - 1, 1);
  }

  private tofAt(j: number) {
    const g = this.grid!;
    return g.tofMin + ((g.tofMax - g.tofMin) * j) / Math.max(g.nTof - 1, 1);
  }

  /** Recomputes the metric and its colour image. */
  private update() {
    const g = this.grid;
    if (!g) return;
    const { of } = METRICS[this.metric];
    this.values = g.c3.map((c3, k) => of(c3, g.vinfArr[k]));
    let min = Infinity;
    this.best = -1;
    this.values.forEach((v, k) => {
      if (v < min) {
        min = v;
        this.best = k;
      }
    });
    // Log colour scale from the minimum to eight times it; anything dearer is the darkest colour.
    const scale = (v: number) => Math.log(v / min) / Math.log(8);

    this.image.width = g.nDep;
    this.image.height = g.nTof;
    const ictx = this.image.getContext("2d")!;
    const img = ictx.createImageData(g.nDep, g.nTof);
    for (let j = 0; j < g.nTof; j++) {
      for (let i = 0; i < g.nDep; i++) {
        const v = this.values[j * g.nDep + i];
        const o = ((g.nTof - 1 - j) * g.nDep + i) * 4; // flight time increases upwards
        if (Number.isFinite(v)) {
          const [r, gr, b] = colour(min > 0 ? scale(v) : 0);
          img.data.set([r, gr, b, 255], o);
        }
      }
    }
    ictx.putImageData(img, 0, 0);

    const levels = min > 0 ? niceLevels(min, 8 * min) : [];
    this.lines = contours(this.values, g.nDep, g.nTof, levels).map((segs, k) => ({ level: levels[k], segs }));
    this.draw();
  }

  private plotRect() {
    const { width, height } = this.canvas.getBoundingClientRect();
    return { x: M.left, y: M.top, w: width - M.left - M.right, h: height - M.top - M.bottom, width, height };
  }

  private cellAt(px: number, py: number): [number, number] | null {
    const g = this.grid;
    if (!g) return null;
    const r = this.plotRect();
    const fx = (px - r.x) / r.w;
    const fy = 1 - (py - r.y) / r.h;
    if (fx < 0 || fx > 1 || fy < 0 || fy > 1) return null;
    return [g.depStart + fx * (g.depEnd - g.depStart), g.tofMin + fy * (g.tofMax - g.tofMin)];
  }

  private hover(px: number, py: number) {
    const cell = this.cellAt(px, py);
    this.canvas.style.cursor = cell ? "crosshair" : "default";
    this.showCell(cell ?? this.picked);
  }

  /** Writes the values at the grid cell nearest a departure date and flight time. */
  private showCell(cell: [number, number] | null) {
    const g = this.grid;
    if (!g || !cell) {
      this.readout.textContent = "";
      return;
    }
    const [dep, tof] = cell;
    const i = Math.round(((dep - g.depStart) / (g.depEnd - g.depStart || 1)) * (g.nDep - 1));
    const j = Math.round(((tof - g.tofMin) / (g.tofMax - g.tofMin || 1)) * (g.nTof - 1));
    const k = j * g.nDep + i;
    const c3 = g.c3[k];
    const va = g.vinfArr[k];
    this.readout.textContent = Number.isFinite(c3)
      ? `Depart ${fmtDate(dep)} · ${tof.toFixed(0)} d · arrive ${fmtDate(dep + tof)}\n` +
        `C3 ${c3.toFixed(1)} km²/s² · v∞ arr ${va.toFixed(2)} km/s`
      : `Depart ${fmtDate(dep)} · ${tof.toFixed(0)} d · no solution`;
  }

  private draw() {
    const g = this.grid;
    const r = this.plotRect();
    const dpr = window.devicePixelRatio || 1;
    this.canvas.width = Math.round(r.width * dpr);
    this.canvas.height = Math.round(r.height * dpr);
    const ctx = this.ctx;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.fillStyle = BG;
    ctx.fillRect(0, 0, r.width, r.height);
    if (!g) return;

    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(this.image, r.x, r.y, r.w, r.h);

    const toX = (dep: number) => r.x + ((dep - g.depStart) / (g.depEnd - g.depStart || 1)) * r.w;
    const toY = (tof: number) => r.y + r.h - ((tof - g.tofMin) / (g.tofMax - g.tofMin || 1)) * r.h;

    // Contours, each labelled once near the middle of its segments.
    const gx = (i: number) => r.x + (i / Math.max(g.nDep - 1, 1)) * r.w;
    const gy = (j: number) => r.y + r.h - (j / Math.max(g.nTof - 1, 1)) * r.h;
    ctx.save();
    ctx.beginPath();
    ctx.rect(r.x, r.y, r.w, r.h);
    ctx.clip();
    ctx.strokeStyle = "#ffffffaa";
    ctx.lineWidth = 0.8;
    ctx.beginPath();
    for (const { segs } of this.lines) {
      for (let k = 0; k + 3 < segs.length; k += 4) {
        ctx.moveTo(gx(segs[k]), gy(segs[k + 1]));
        ctx.lineTo(gx(segs[k + 2]), gy(segs[k + 3]));
      }
    }
    ctx.stroke();
    ctx.font = "10px system-ui, sans-serif";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    for (const { level, segs } of this.lines) {
      if (segs.length < 4) continue;
      // Label where the contour crosses a diagonal ray up and right from the minimum, so nested
      // rings label in a line without overlapping.
      const bi = this.best % g.nDep;
      const bj = Math.floor(this.best / g.nDep);
      let k = 0;
      let off = Infinity;
      for (let s = 0; s + 3 < segs.length; s += 4) {
        const dx = ((segs[s] + segs[s + 2]) / 2 - bi) / g.nDep;
        const dy = ((segs[s + 1] + segs[s + 3]) / 2 - bj) / g.nTof;
        const d = Math.abs(Math.atan2(dy, dx) - Math.PI / 4);
        if (d < off) {
          off = d;
          k = s;
        }
      }
      const x = gx((segs[k] + segs[k + 2]) / 2);
      const y = gy((segs[k + 1] + segs[k + 3]) / 2);
      const text = level >= 10 ? level.toFixed(0) : level.toFixed(1);
      const w = ctx.measureText(text).width + 4;
      ctx.fillStyle = BG + "cc";
      ctx.fillRect(x - w / 2, y - 6, w, 12);
      ctx.fillStyle = "#ffffff";
      ctx.fillText(text, x, y);
    }
    ctx.restore();

    // Axes: months along x, flight time up y.
    ctx.font = FONT;
    ctx.fillStyle = INK;
    ctx.strokeStyle = INK + "60";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    const start = dateOf(g.depStart);
    const months = Math.ceil((g.depEnd - g.depStart) / 30.44);
    const every = Math.max(1, Math.ceil(months / 6));
    for (let m = 1; m <= months; m++) {
      const d = new Date(Date.UTC(start.getUTCFullYear(), start.getUTCMonth() + m, 1));
      const x = toX((d.getTime() - Date.UTC(2000, 0, 1, 12)) / 86_400_000);
      if (x > r.x + r.w) break;
      ctx.beginPath();
      ctx.moveTo(x, r.y + r.h);
      ctx.lineTo(x, r.y + r.h + 3);
      ctx.stroke();
      if (m % every === 0)
        ctx.fillText(d.toLocaleDateString("en-GB", { month: "short", year: "2-digit", timeZone: "UTC" }), x, r.y + r.h + 5);
    }
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    const range = g.tofMax - g.tofMin;
    const tick = [10, 20, 50, 100, 200, 500, 1000, 2000].find((s) => range / s <= 6) ?? 5000;
    for (let t = Math.ceil(g.tofMin / tick) * tick; t <= g.tofMax; t += tick) {
      const y = toY(t);
      ctx.beginPath();
      ctx.moveTo(r.x - 3, y);
      ctx.lineTo(r.x, y);
      ctx.stroke();
      ctx.fillText(`${t} d`, r.x - 5, y);
    }

    // Cheapest cell, and the picked one.
    if (this.best >= 0) {
      const x = toX(this.depAt(this.best % g.nDep));
      const y = toY(this.tofAt(Math.floor(this.best / g.nDep)));
      ctx.beginPath();
      ctx.arc(x, y, 5, 0, Math.PI * 2);
      ctx.strokeStyle = "#ffffff";
      ctx.lineWidth = 1.5;
      ctx.stroke();
    }
    if (this.picked) {
      const x = toX(this.picked[0]);
      const y = toY(this.picked[1]);
      ctx.strokeStyle = "#ff5d8f";
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(x, r.y);
      ctx.lineTo(x, r.y + r.h);
      ctx.moveTo(r.x, y);
      ctx.lineTo(r.x + r.w, y);
      ctx.stroke();
    }
    ctx.lineWidth = 1;
  }
}
