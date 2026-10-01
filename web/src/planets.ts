import type { Planet } from "./engine/engine.js";
import type { Item } from "./plot/items.js";

export const COLOURS: Record<string, string> = {
  Mercury: "#a8a29e",
  Venus: "#e8c07a",
  Earth: "#4a90d9",
  Mars: "#d0603a",
  Jupiter: "#d9a066",
  Saturn: "#e3cf8f",
  Uranus: "#8fd3dc",
  Neptune: "#4f6fd6",
};

/** Orbit lines for each planet. */
export function orbitLayer(bodies: Planet[]): Item[] {
  return bodies.map((p) => ({ kind: "line", xy: p.orbit_xy(360), color: COLOURS[p.name] + "80" }));
}

/** Planet markers at `day` (days after J2000). */
export function planetMarkers(bodies: Planet[], day: number): Item[] {
  return bodies.map((p) => {
    const [x, y] = p.position_at(day);
    return {
      kind: "marker",
      id: p.name,
      x,
      y,
      color: COLOURS[p.name],
      size: p.a_au > 4 ? 7 : 5,
      label: p.name,
    };
  });
}

/** Fills `info` with a planet's orbital elements. */
export function showElements(info: HTMLElement, p: Planet) {
  const deg = (v: number) => `${v.toFixed(3)}°`;
  const rows: [string, string][] = [
    ["Semi-major axis a", `${p.a_au.toFixed(4)} AU`],
    ["Eccentricity e", p.e.toFixed(5)],
    ["Inclination i", deg(p.i_deg)],
    ["Ascending node Ω", deg(p.raan_deg)],
    ["Longitude of perihelion ϖ", deg(p.lon_peri_deg)],
    ["Mean longitude L", deg(p.mean_lon_deg)],
    ["Period", `${p.period_days.toFixed(1)} d (${(p.period_days / 365.25).toFixed(2)} yr)`],
  ];
  info.innerHTML = `<h2>${p.name}</h2><dl>${rows
    .map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`)
    .join("")}</dl><p class="hint">J2000 elements, ecliptic frame.</p>`;
}
