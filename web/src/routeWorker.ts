/** Runs the flyby-sequence search off the main thread. */
import init, { find_routes } from "./engine/engine.js";
import type { FoundRoute, SearchRequest, WorkerMessage } from "./routeFinder.js";

const ready = init();
const send = (msg: WorkerMessage) => postMessage(msg);

onmessage = async (e: MessageEvent<SearchRequest>) => {
  await ready;
  const q = e.data;
  const started = performance.now();
  const found = find_routes(
    q.from, q.to, q.launchStart, q.launchEnd, q.maxFlybys, q.capture,
    q.timeCost, q.maxTof, q.step, q.margin, new Uint32Array(q.flybys),
    (status: string) => send({ type: "progress", status }),
  );
  const routes: FoundRoute[] = found.map((r) => {
    const out = {
      bodies: Array.from(r.bodies), dates: Array.from(r.dates), c3: r.c3, departure: r.departure,
      flybys: Array.from(r.flybys), arrival: r.arrival, dv: r.dv, score: r.score,
    };
    r.free();
    return out;
  });
  send({ type: "done", routes, ms: performance.now() - started });
};
