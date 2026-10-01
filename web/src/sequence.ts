/** Planet letters for compact sequences like EVVEJS. M is Mars; Mercury has to be written Me. */
const LETTERS: Record<string, string> = { E: "Earth", V: "Venus", M: "Mars", J: "Jupiter", S: "Saturn", U: "Uranus", N: "Neptune" };

/**
 * Parses the planets to fly past from text such as "VVEJ", "EVVEJS", "Venus Venus Earth Jupiter"
 * or "E, V, Me". A word is read as a planet name if it starts one (at least two letters, so "Me"
 * is Mercury); otherwise as one letter per planet. If the text ends with the target, it is taken
 * as the whole route and the departure and target are dropped.
 *
 * Returns planet indices (empty for "search all"), or an error message.
 */
export function parseFlybys(text: string, names: string[], from: number, to: number): number[] | string {
  const words = text.trim().split(/[\s,>–-]+/).filter(Boolean);
  const seq: number[] = [];
  for (const word of words) {
    const w = word.toLowerCase();
    const named = w.length >= 2 ? names.findIndex((n) => n.toLowerCase().startsWith(w)) : -1;
    if (named >= 0) {
      seq.push(named);
      continue;
    }
    for (const ch of word.toUpperCase()) {
      const k = names.indexOf(LETTERS[ch]);
      if (k < 0) return `Unknown planet "${ch}" in "${word}".`;
      seq.push(k);
    }
  }
  if (seq.length && seq[seq.length - 1] === to) {
    seq.pop();
    if (seq[0] === from) seq.shift();
  }
  if (seq.includes(to)) return `${names[to]} is the target, so it can't be flown past.`;
  return seq;
}
