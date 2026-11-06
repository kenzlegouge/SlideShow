"""Local web server of the viewer (Python standard library only).

Serves the viewer (package folder app/) at /, the bundle at /bundle/, and stores the labels typed in the page:
    <labels>/cluster_labels.csv    cluster, label, note, n_cells   <- slideshow.read_labels / apply_labels
    <labels>/cluster_labels.json   same, plus the clustering fingerprint the page checks on reload
    <labels>/history/              copy of the CSV as it was when this server session first overwrote it
Listens on 127.0.0.1 only.
"""
import csv
import json
import os
import shutil
import threading
import time
import webbrowser
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from .labels import labels_dir

APP = Path(__file__).resolve().parent / "app"


def serve(bundle, port=8765, labels=None, open_browser=True, block=True):
    """Serve the viewer for `bundle` at http://127.0.0.1:<port>/.

    labels        folder for the labels (default <bundle>/labels)
    block=False   run in a background thread (e.g. from a notebook); returns the server, stop it with .shutdown()
    """
    bundle = Path(bundle).resolve()
    meta_file = bundle / "meta.json"
    if not meta_file.exists():
        raise FileNotFoundError(f"{bundle} has no meta.json: build it with slideshow.build first")
    meta = json.loads(meta_file.read_text(encoding="utf-8"))
    if meta.get("format") != "slideshow-bundle":
        raise ValueError(f"{bundle} is not a SlideShow bundle")
    ann = labels_dir(bundle, labels).resolve()
    lock = threading.Lock()
    state = {"backed_up": False}

    class Handler(SimpleHTTPRequestHandler):
        def translate_path(self, path):
            rel = path.split("?", 1)[0].split("#", 1)[0].lstrip("/")
            if rel.startswith("bundle/"):
                root, rel = bundle, rel[len("bundle/"):]
            else:
                root, rel = APP, rel or "index.html"
            target = (root / rel).resolve()
            return str(target) if target.is_relative_to(root) else str(root / "__forbidden__")

        def end_headers(self):
            self.send_header("Cache-Control", "no-cache")      # always revalidate: the bundle can be rebuilt
            super().end_headers()

        def send_json(self, obj, status=200):
            body = json.dumps(obj).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            if self.path.split("?", 1)[0] == "/api/annotations":
                f = ann / "cluster_labels.json"
                return self.send_json(json.loads(f.read_text(encoding="utf-8")) if f.exists() else {})
            super().do_GET()

        def do_POST(self):
            if self.path != "/api/annotations":
                return self.send_json({"error": "not found"}, 404)
            try:
                data = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))))
                rows = data["rows"]
                if not isinstance(rows, list) or not all(isinstance(r, dict) and "cluster" in r for r in rows):
                    raise ValueError("rows must be a list of {cluster, label, note, n_cells}")
            except (ValueError, KeyError) as e:
                return self.send_json({"error": str(e)}, 400)
            with lock:
                ann.mkdir(parents=True, exist_ok=True)
                csv_path = ann / "cluster_labels.csv"
                if csv_path.exists() and not state["backed_up"]:
                    (ann / "history").mkdir(exist_ok=True)
                    shutil.copy2(csv_path, ann / "history" / f"cluster_labels_{time.strftime('%Y%m%d_%H%M%S')}.csv")
                state["backed_up"] = True
                tmp = csv_path.with_suffix(".csv.tmp")
                with open(tmp, "w", newline="", encoding="utf-8") as fh:
                    w = csv.writer(fh)
                    w.writerow(["cluster", "label", "note", "n_cells"])
                    for r in rows:
                        w.writerow([r["cluster"], r.get("label", ""), r.get("note", ""), r.get("n_cells", "")])
                os.replace(tmp, csv_path)
                tmp = ann / "cluster_labels.json.tmp"
                tmp.write_text(json.dumps(data, indent=1), encoding="utf-8")
                os.replace(tmp, ann / "cluster_labels.json")
            self.send_json({"ok": True, "csv": str(csv_path)})

        def log_message(self, fmt, *a):
            if self.command == "POST" or (a and str(a[1])[:1] in "45"):
                super().log_message(fmt, *a)

    url = f"http://127.0.0.1:{port}/"
    srv = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    print(f"SlideShow: {url}\n  bundle {bundle}\n  labels {ann}")
    if open_browser:
        threading.Timer(0.5, webbrowser.open, [url]).start()
    if not block:
        threading.Thread(target=srv.serve_forever, daemon=True).start()
        return srv
    print("Ctrl+C to stop.")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.server_close()
