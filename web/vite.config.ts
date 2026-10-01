import { defineConfig } from "vite";

// ES2022 for top-level await in main.ts.
export default defineConfig({
  build: { target: "es2022" },
});
