"""Static server for the browser tests: cross-origin isolation headers (SharedArrayBuffer), wasm MIME type,
and POST /result (the page posts its log; saved under web/results/). Port from $PORT or argv[1]."""
import http.server
import os
import sys
import time

os.chdir(os.path.dirname(os.path.abspath(__file__)))


class Handler(http.server.SimpleHTTPRequestHandler):
    extensions_map = {**http.server.SimpleHTTPRequestHandler.extensions_map, ".wasm": "application/wasm"}

    def end_headers(self):
        self.send_header("Cross-Origin-Opener-Policy", "same-origin")
        self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def do_POST(self):
        body = self.rfile.read(int(self.headers.get("Content-Length", 0)))
        os.makedirs("results", exist_ok=True)
        browser = "firefox" if b"Firefox" in body[:400] else "chromium"
        with open(f"results/{browser}-{int(time.time())}.txt", "wb") as f:
            f.write(body)
        self.send_response(204)
        self.end_headers()


port = int(os.environ.get("PORT", sys.argv[1] if len(sys.argv) > 1 else 8766))
http.server.ThreadingHTTPServer(("127.0.0.1", port), Handler).serve_forever()
