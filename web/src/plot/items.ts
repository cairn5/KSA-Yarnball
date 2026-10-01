/** Things the plot can draw. Coordinates are heliocentric ecliptic x, y in AU. */

export type Shape = "circle" | "triangle" | "square" | "star" | "cross";

/** A polyline, e.g. an orbit or a trajectory leg. `xy` is interleaved x, y. */
export interface Line {
  kind: "line";
  xy: ArrayLike<number>;
  color: string;
  width?: number;
  dash?: number[];
}

/** A point marker, e.g. a planet or an event. Markers with an `id` can be clicked. */
export interface Marker {
  kind: "marker";
  x: number;
  y: number;
  color: string;
  shape?: Shape;
  size?: number;
  id?: string;
  /** Text under the marker. */
  label?: string;
  /** Short text in a circle beside the marker, e.g. an event number. */
  tag?: string;
}

export type Item = Line | Marker;

/** Maps AU to canvas pixels. */
export type Project = (x: number, y: number) => [number, number];

export const INK = "#c9d3e6";
export const BG = "#0b1020";
export const FONT = "11px system-ui, sans-serif";

export const markerSize = (m: Marker) => m.size ?? 5;

export function drawItem(ctx: CanvasRenderingContext2D, item: Item, project: Project) {
  if (item.kind === "line") drawLine(ctx, item, project);
  else drawMarker(ctx, item, project);
}

function drawLine(ctx: CanvasRenderingContext2D, line: Line, project: Project) {
  const { xy } = line;
  ctx.beginPath();
  for (let i = 0; i + 1 < xy.length; i += 2) {
    const [sx, sy] = project(xy[i], xy[i + 1]);
    if (i === 0) ctx.moveTo(sx, sy);
    else ctx.lineTo(sx, sy);
  }
  ctx.strokeStyle = line.color;
  ctx.lineWidth = line.width ?? 1;
  ctx.setLineDash(line.dash ?? []);
  ctx.stroke();
  ctx.setLineDash([]);
}

function drawMarker(ctx: CanvasRenderingContext2D, m: Marker, project: Project) {
  const [x, y] = project(m.x, m.y);
  const s = markerSize(m);
  ctx.beginPath();
  switch (m.shape ?? "circle") {
    case "circle":
      ctx.arc(x, y, s, 0, Math.PI * 2);
      break;
    case "square":
      ctx.rect(x - s, y - s, 2 * s, 2 * s);
      break;
    case "triangle":
      ctx.moveTo(x, y - 1.2 * s);
      ctx.lineTo(x + 1.1 * s, y + 0.8 * s);
      ctx.lineTo(x - 1.1 * s, y + 0.8 * s);
      ctx.closePath();
      break;
    case "star":
      for (let k = 0; k < 10; k++) {
        const r = k % 2 === 0 ? 1.5 * s : 0.6 * s;
        const a = -Math.PI / 2 + (k * Math.PI) / 5;
        ctx.lineTo(x + r * Math.cos(a), y + r * Math.sin(a));
      }
      ctx.closePath();
      break;
    case "cross":
      ctx.moveTo(x - s, y - s);
      ctx.lineTo(x + s, y + s);
      ctx.moveTo(x + s, y - s);
      ctx.lineTo(x - s, y + s);
      ctx.strokeStyle = m.color;
      ctx.lineWidth = 1.5;
      ctx.stroke();
      break;
  }
  if (m.shape !== "cross") {
    ctx.fillStyle = m.color;
    ctx.fill();
    ctx.strokeStyle = BG;
    ctx.lineWidth = 1;
    ctx.stroke();
  }

  ctx.font = FONT;
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";
  if (m.label) {
    ctx.fillStyle = INK;
    ctx.fillText(m.label, x, y + s + 10);
  }
  if (m.tag) {
    const tx = x + s + 10;
    const ty = y - s - 10;
    ctx.beginPath();
    ctx.arc(tx, ty, 9, 0, Math.PI * 2);
    ctx.fillStyle = BG;
    ctx.fill();
    ctx.strokeStyle = INK;
    ctx.stroke();
    ctx.fillStyle = INK;
    ctx.fillText(m.tag, tx, ty + 0.5);
  }
}
