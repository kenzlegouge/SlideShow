"""Browser smoke test: serve the synthetic bundle, open it in headless Chrome (software WebGL) over the DevTools
protocol, and check the viewer loads, colours, computes markers, filters and draws polygons without a JS error.
Skipped when Chrome or aiohttp is not available. Screenshots go to the pytest tmp folder (printed with -s)."""
import asyncio
import base64
import json
import os
import shutil
import socket
import subprocess
import time
from pathlib import Path

import pytest

aiohttp = pytest.importorskip("aiohttp")
CHROME = os.environ.get("SLIDESHOW_CHROME") or next(
    (p for p in [r"C:\Program Files\Google\Chrome\Application\chrome.exe", shutil.which("google-chrome"),
                 shutil.which("chromium"), shutil.which("chrome")] if p and Path(p).exists()), None)
pytestmark = pytest.mark.skipif(CHROME is None, reason="Chrome not found (set SLIDESHOW_CHROME)")


def free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def test_viewer(bundle, tmp_path):
    from slideshow import serve
    port, debug_port = free_port(), free_port()
    srv = serve(bundle, port=port, open_browser=False, block=False)
    proc = subprocess.Popen([CHROME, "--headless=new", f"--remote-debugging-port={debug_port}", f"--user-data-dir={tmp_path / 'profile'}",
                             "--window-size=1600,900", "--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--no-first-run", "about:blank"],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        errors, report = asyncio.run(drive(f"http://127.0.0.1:{port}/", debug_port, tmp_path))
    finally:
        proc.terminate()
        srv.shutdown()
        srv.server_close()
    print("\nviewer:", json.dumps(report, indent=1), "\nscreenshots:", tmp_path)
    assert not errors, errors
    assert report["subtitle"].startswith("synthetic · 4,800 cells · 40 genes · 5 clusters")
    assert report["top_marker"] in [f"G{j}" for j in range(8)]          # cluster 0 expresses G0..G7
    assert report["replicates_evaluated"] == "6 cores evaluated"
    assert report["kept_after_filter"] < 4800 and report["kept_after_filter"] > 1000
    assert report["polygons_drawn"]


async def drive(url, debug_port, shots):
    errors, pending, n = [], {}, [0]
    async with aiohttp.ClientSession() as http:
        for _ in range(100):
            try:
                tabs = await (await http.get(f"http://127.0.0.1:{debug_port}/json")).json()
                break
            except Exception:
                await asyncio.sleep(0.2)
        page = next(t for t in tabs if t["type"] == "page")
        async with http.ws_connect(page["webSocketDebuggerUrl"], max_msg_size=0) as ws:

            async def reader():
                async for msg in ws:
                    d = json.loads(msg.data)
                    if "id" in d and d["id"] in pending:
                        pending.pop(d["id"]).set_result(d)
                    elif d.get("method") == "Runtime.exceptionThrown":
                        ex = d["params"]["exceptionDetails"]
                        errors.append(ex.get("exception", {}).get("description") or ex["text"])
                    elif d.get("method") == "Runtime.consoleAPICalled" and d["params"]["type"] == "error":
                        errors.append("console.error: " + " ".join(str(a.get("value", a.get("description"))) for a in d["params"]["args"]))

            async def send(method, **params):
                n[0] += 1
                fut = asyncio.get_running_loop().create_future()
                pending[n[0]] = fut
                await ws.send_str(json.dumps({"id": n[0], "method": method, "params": params}))
                r = await asyncio.wait_for(fut, 60)
                if "error" in r:
                    raise RuntimeError(r["error"])
                return r["result"]

            async def js(expr):
                r = await send("Runtime.evaluate", expression=expr, awaitPromise=True, returnByValue=True)
                if "exceptionDetails" in r:
                    raise RuntimeError(r["exceptionDetails"])
                return r["result"].get("value")

            async def wait(expr, timeout=60):
                t0 = time.time()
                while time.time() - t0 < timeout:
                    if await js(expr):
                        return
                    await asyncio.sleep(0.25)
                raise TimeoutError(f"{expr}; errors: {errors}")

            async def shot(name):
                await asyncio.sleep(0.5)
                r = await send("Page.captureScreenshot", format="png")
                (shots / f"{name}.png").write_bytes(base64.b64decode(r["data"]))

            asyncio.ensure_future(reader())
            await send("Runtime.enable")
            await send("Page.enable")
            await send("Emulation.setDeviceMetricsOverride", width=1600, height=900, deviceScaleFactor=1, mobile=False)
            await send("Page.navigate", url=url)
            await wait("document.getElementById('loading') === null || /Could not/.test(document.getElementById('loading').innerText)")
            loading = await js("document.getElementById('loading') && document.getElementById('loading').innerText")
            if loading:
                return [loading], {}
            await wait("LAYOUT.tissueReady")
            report = {"subtitle": await js("document.getElementById('subtitle').textContent")}
            await shot("01_loaded")
            # select cluster 0: marker table
            await js("select(0), true")
            await asyncio.sleep(1)
            report["top_marker"] = await js("document.querySelector('#markers tbody tr td').textContent")
            report["replicates_evaluated"] = await js("document.getElementById('markertitle').textContent.split(' · ')[1]")
            await shot("02_cluster0")
            # colour by gene, by the summary category, by a continuous covariate
            for v in ["cat:auto_type", "cat:region", "cont:score", "label"]:
                await js(f"(() => {{ const s = document.getElementById('colorby'); s.value = '{v}'; s.dispatchEvent(new Event('change')); return true; }})()")
                await asyncio.sleep(0.4)
            await shot("03_colour_score")
            # filter: hide 'Tumor'
            await js("addFilter('region')")
            await asyncio.sleep(0.5)
            await js("document.querySelector('#filters .flevels input[type=checkbox]:nth-of-type(1)') ? true : true")
            await js("(() => { const box = [...document.querySelectorAll('#filters .flevels label')].find(l => l.textContent.includes('Tumor')).querySelector('input'); box.click(); return true; })()")
            await asyncio.sleep(1)
            report["kept_after_filter"] = await js("S.keep ? S.keep.reduce((a, b) => a + b, 0) : -1")
            report["filter_count"] = await js("document.getElementById('filtercount').textContent")
            # tissue view: open a sample and zoom into polygons
            await js("LAYOUT.umap = false; LAYOUT.side = false; applyLayout(); true")
            await js("SP.openSample('S2_core1'), true")
            for _ in range(6):
                await send("Input.dispatchMouseEvent", type="mouseWheel", x=800, y=500, deltaX=0, deltaY=-400)
                await asyncio.sleep(0.3)
            await wait("!document.getElementById('spinfo').innerText.includes('loading')")
            await asyncio.sleep(1)
            report["tissue_info"] = await js("document.getElementById('spinfo').innerText")
            report["polygons_drawn"] = not any(w in report["tissue_info"] for w in ("zoom in", "no segmentation"))
            await shot("04_tissue_polygons")
            return errors, report
