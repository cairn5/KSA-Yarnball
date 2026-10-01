/**
 * Contour lines by marching squares. `values` is row-major (`j * nx + i`); NaN cells are skipped.
 * Returns, per level, line segments as interleaved grid coordinates [i1, j1, i2, j2, ...].
 */
export function contours(values: ArrayLike<number>, nx: number, ny: number, levels: number[]): number[][] {
  return levels.map((level) => {
    const segs: number[] = [];
    const pts: number[] = [];
    for (let j = 0; j + 1 < ny; j++) {
      for (let i = 0; i + 1 < nx; i++) {
        // Corners anticlockwise from (i, j).
        const c = [
          [i, j, values[j * nx + i]],
          [i + 1, j, values[j * nx + i + 1]],
          [i + 1, j + 1, values[(j + 1) * nx + i + 1]],
          [i, j + 1, values[(j + 1) * nx + i]],
        ];
        if (!c.every((p) => Number.isFinite(p[2]))) continue;
        pts.length = 0;
        for (let k = 0; k < 4; k++) {
          const [xa, ya, va] = c[k];
          const [xb, yb, vb] = c[(k + 1) % 4];
          if (va < level !== vb < level) {
            const t = (level - va) / (vb - va);
            pts.push(xa + t * (xb - xa), ya + t * (yb - ya));
          }
        }
        // Two crossings make one segment; four (a saddle) make two.
        segs.push(...pts);
      }
    }
    return segs;
  });
}

/** Round numbers between `lo` and `hi`, about `target` of them. */
export function niceLevels(lo: number, hi: number, target = 7): number[] {
  const mantissas = [1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8];
  const all: number[] = [];
  for (let e = Math.floor(Math.log10(lo)) - 1; e <= Math.ceil(Math.log10(hi)); e++) {
    for (const m of mantissas) {
      const v = m * 10 ** e;
      if (v > lo && v <= hi) all.push(v);
    }
  }
  const every = Math.max(1, Math.ceil(all.length / target));
  return all.filter((_, k) => k % every === 0);
}
