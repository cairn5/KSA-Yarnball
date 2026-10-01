/** Dates as days after J2000 (2000-01-01 12:00), the engine's time scale. UTC vs TT (~1 min) is ignored. */

const J2000_MS = Date.UTC(2000, 0, 1, 12);
const DAY_MS = 86_400_000;

export const daysOf = (d: Date) => (d.getTime() - J2000_MS) / DAY_MS;
export const dateOf = (days: number) => new Date(J2000_MS + days * DAY_MS);

/** YYYY-MM-DD. */
export const fmtDate = (days: number) => dateOf(days).toISOString().slice(0, 10);

/** Days after J2000 for a YYYY-MM-DD string, at 00:00 UTC. */
export const parseDate = (s: string) => daysOf(new Date(`${s}T00:00:00Z`));
