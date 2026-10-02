import { BG, FONT, INK, drawItem, markerSize, type Item, type Marker } from "./items.js";

export interface PlotOptions {
  /** Radial scale length, AU: radius is drawn as ln(1 + r / r0). */
  r0?: number;
  /** Radii of the dotted reference rings, AU. */
  gridAu?: number[];
  /** Radius that fits the screen at the default zoom, AU. */
  fitAu?: number;
}

/**
 * Plan view of the ecliptic from the north, on a canvas. Radius is log-scaled (not to scale),
 * directions are true. Drag to pan, wheel or pinch to zoom, double-click to reset.
 */
export class Plot {
  /** Id of the selected marker, ringed when drawn. */
  selected: string | null = null;
  /** Called when the user clicks a marker with an id, or empty space (null). */
  onSelect: (id: string | null) => void = () => {};

  private readonly ctx: CanvasRenderingContext2D;
  private readonly layers = new Map<string, Item[]>();
  private readonly r0: number;
  private readonly gridAu: number[];
  private readonly fitAu: number;
  private width = 0;
  private height = 0;
  private k = 0; // pixels per unit of ln(1 + r / r0)
  private panX = 0;
  private panY = 0;
  private queued = false;

  constructor(private readonly canvas: HTMLCanvasElement, opts: PlotOptions = {}) {
    this.ctx = canvas.getContext("2d")!;
    this.r0 = opts.r0 ?? 0.3;
    this.gridAu = opts.gridAu ?? [0.5, 1, 2, 5, 10, 20, 30];
    this.fitAu = opts.fitAu ?? 32;
    new ResizeObserver(() => this.resize()).observe(canvas);
    this.resize();
    this.bindInput();
  }

  /** Replaces a named layer. Layers draw in the order they were first set. */
  setLayer(name: string, items: Item[]) {
    this.layers.set(name, items);
    this.render();
  }

  select(id: string | null) {
    this.selected = id;
    this.onSelect(id);
    this.render();
  }

  /** Heliocentric ecliptic x, y in AU to canvas pixels. */
  project = (x: number, y: number): [number, number] => {
    const r = Math.hypot(x, y);
    const s = r === 0 ? 0 : (Math.log1p(r / this.r0) / r) * this.k;
    return [this.width / 2 + this.panX + x * s, this.height / 2 + this.panY - y * s];
  };

  /** Back to the default view. */
  reset() {
    this.k = (0.46 * Math.min(this.width, this.height)) / Math.log1p(this.fitAu / this.r0);
    this.panX = this.panY = 0;
    this.render();
  }

  render() {
    if (this.queued) return;
    this.queued = true;
    requestAnimationFrame(() => {
      this.queued = false;
      this.draw();
    });
  }

  private resize() {
    const dpr = window.devicePixelRatio || 1;
    const { width, height } = this.canvas.getBoundingClientRect();
    this.width = width;
    this.height = height;
    this.canvas.width = Math.round(width * dpr);
    this.canvas.height = Math.round(height * dpr);
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    if (this.k === 0) this.reset();
    else this.render();
  }

  private draw() {
    const { ctx } = this;
    ctx.fillStyle = BG;
    ctx.fillRect(0, 0, this.width, this.height);
    this.drawWatermark();

    const [cx, cy] = this.project(0, 0);
    ctx.font = FONT;
    ctx.textAlign = "center";
    ctx.textBaseline = "bottom";
    ctx.setLineDash([1, 4]);
    ctx.strokeStyle = "#ffffff30";
    ctx.lineWidth = 1;
    for (const au of this.gridAu) {
      const r = Math.log1p(au / this.r0) * this.k;
      ctx.beginPath();
      ctx.arc(cx, cy, r, 0, Math.PI * 2);
      ctx.stroke();
      ctx.fillStyle = "#ffffff60";
      ctx.fillText(`${au} AU`, cx, cy - r - 2);
    }
    ctx.setLineDash([]);

    ctx.beginPath();
    ctx.arc(cx, cy, 8, 0, Math.PI * 2);
    ctx.fillStyle = "#f5c542";
    ctx.fill();

    for (const items of this.layers.values()) for (const item of items) drawItem(ctx, item, this.project);

    const sel = this.selected && [...this.markers()].find((m) => m.id === this.selected);
    if (sel) {
      const [x, y] = this.project(sel.x, sel.y);
      ctx.beginPath();
      ctx.arc(x, y, markerSize(sel) + 6, 0, Math.PI * 2);
      ctx.strokeStyle = "#ffd166";
      ctx.lineWidth = 2;
      ctx.stroke();
    }

    ctx.fillStyle = INK + "99";
    ctx.textAlign = "left";
    ctx.textBaseline = "bottom";
    ctx.fillText(`radius ln(1 + r / ${this.r0} AU): not to scale, directions true`, 12, this.height - 10);
  }

  /** Faint widely-spaced title fixed to the screen centre, behind everything. */
  private drawWatermark() {
    const { ctx } = this;
    const size = Math.min(this.width / 8, this.height / 5);
    const cx = this.width / 2;
    const cy = this.height / 2;
    ctx.save();
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillStyle = "#ffffff0d";
    ctx.font = `200 ${size}px system-ui, sans-serif`;
    ctx.letterSpacing = `${size * 0.35}px`;
    // Trailing letter spacing pushes the text left; shift it back by half.
    ctx.fillText("YARNBALL", cx + size * 0.175, cy);
    const small = size * 0.16;
    ctx.font = `300 ${small}px system-ui, sans-serif`;
    ctx.letterSpacing = `${small * 0.6}px`;
    ctx.fillText("DISTANCES LOG SCALE", cx + small * 0.3, cy + size * 0.75);
    ctx.restore();
  }

  private *markers(): Generator<Marker> {
    for (const items of this.layers.values())
      for (const item of items) if (item.kind === "marker") yield item;
  }

  /** Topmost clickable marker within reach of a canvas point. */
  private hit(px: number, py: number): string | null {
    let best: string | null = null;
    let bestD = Infinity;
    for (const m of this.markers()) {
      if (!m.id) continue;
      const [x, y] = this.project(m.x, m.y);
      const d = Math.hypot(x - px, y - py);
      if (d < markerSize(m) + 8 && d <= bestD) {
        best = m.id;
        bestD = d;
      }
    }
    return best;
  }

  private zoomAt(px: number, py: number, factor: number) {
    const k = Math.min(Math.max(this.k * factor, 20), 50_000);
    const f = k / this.k;
    // Keep the point under the cursor fixed.
    this.panX = px - this.width / 2 - (px - this.width / 2 - this.panX) * f;
    this.panY = py - this.height / 2 - (py - this.height / 2 - this.panY) * f;
    this.k = k;
    this.render();
  }

  private bindInput() {
    const c = this.canvas;
    const pointers = new Map<number, { x: number; y: number }>();
    let moved = 0;

    c.addEventListener("pointerdown", (e) => {
      c.setPointerCapture(e.pointerId);
      pointers.set(e.pointerId, { x: e.offsetX, y: e.offsetY });
      moved = 0;
    });

    c.addEventListener("pointermove", (e) => {
      const p = pointers.get(e.pointerId);
      if (!p) {
        c.style.cursor = this.hit(e.offsetX, e.offsetY) ? "pointer" : "grab";
        return;
      }
      if (pointers.size === 2) {
        // Pinch: zoom by the change in finger spacing, about their midpoint.
        const [a, b] = [...pointers.values()];
        const before = Math.hypot(a.x - b.x, a.y - b.y);
        p.x = e.offsetX;
        p.y = e.offsetY;
        const after = Math.hypot(a.x - b.x, a.y - b.y);
        if (before > 0) this.zoomAt((a.x + b.x) / 2, (a.y + b.y) / 2, after / before);
        moved = Infinity;
        return;
      }
      const dx = e.offsetX - p.x;
      const dy = e.offsetY - p.y;
      p.x = e.offsetX;
      p.y = e.offsetY;
      moved += Math.abs(dx) + Math.abs(dy);
      if (moved > 4) {
        c.style.cursor = "grabbing";
        this.panX += dx;
        this.panY += dy;
        this.render();
      }
    });

    const end = (e: PointerEvent) => {
      if (!pointers.delete(e.pointerId)) return;
      if (moved <= 4 && e.type === "pointerup") this.select(this.hit(e.offsetX, e.offsetY));
      c.style.cursor = "grab";
    };
    c.addEventListener("pointerup", end);
    c.addEventListener("pointercancel", end);

    c.addEventListener(
      "wheel",
      (e) => {
        e.preventDefault();
        this.zoomAt(e.offsetX, e.offsetY, Math.exp(-e.deltaY * 0.0015));
      },
      { passive: false },
    );
    c.addEventListener("dblclick", () => this.reset());
  }
}
