"""Start the loopback-only NeoScorer multi-patient dashboard."""

import argparse
from contextlib import asynccontextmanager
import csv
import io
import json
from pathlib import Path
import socket
import threading
import time
import urllib.request
import webbrowser
import uvicorn
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import FileResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles
from starlette.middleware.trustedhost import TrustedHostMiddleware
from dashboard_backend.jobs import JobManager, RunRequest
from dashboard_backend.reports import ReportStore, filter_rows

BASE = Path(__file__).resolve().parent


def create_app(base=BASE):
    store = ReportStore(base / "reports")
    jobs = JobManager(base, store.directory)

    @asynccontextmanager
    async def lifespan(application):
        yield
        jobs.close()

    app = FastAPI(title="NeoScorer Dashboard", lifespan=lifespan)
    app.state.store, app.state.jobs = store, jobs
    app.add_middleware(
        TrustedHostMiddleware, allowed_hosts=["127.0.0.1", "localhost", "testserver"]
    )

    @app.middleware("http")
    async def local_guard(request: Request, call_next):
        if request.method == "POST":
            origin = request.headers.get("origin")
            if (
                origin
                and origin != f"{request.url.scheme}://{request.headers.get('host')}"
            ):
                return JSONResponse(
                    {"detail": "Cross-origin submissions are not allowed."},
                    status_code=403,
                )
            if request.headers.get("x-neoscorer") != "dashboard":
                return JSONResponse(
                    {"detail": "Missing dashboard request header."}, status_code=403
                )
            body = bytearray()
            async for chunk in request.stream():
                body.extend(chunk)
                if len(body) > 20_000_000:
                    return JSONResponse(
                        {"detail": "Input exceeds the 20 MB limit."}, status_code=413
                    )
            request._body = bytes(body)
        response = await call_next(request)
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["Content-Security-Policy"] = (
            "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; frame-ancestors 'none'; base-uri 'none'; form-action 'self'"
        )
        response.headers["Cache-Control"] = (
            "no-store" if request.url.path.startswith("/api") else "no-cache"
        )
        return response

    app.mount(
        "/static", StaticFiles(directory=BASE / "dashboard_static"), name="static"
    )

    @app.get("/")
    def index():
        return FileResponse(BASE / "dashboard_static" / "index.html")

    @app.get("/api/health")
    def health():
        return {
            "status": "ok",
            "application": "neoscorer-dashboard",
            "patients": len(store.discover()),
        }

    @app.get("/api/patients")
    def patients():
        return store.discover()

    @app.get("/api/patients/{patient}")
    def report(patient: str):
        return store.read(patient)

    @app.get("/api/patients/{patient}/export")
    def export(
        patient: str,
        q: str = "",
        tier: int = Query(0, ge=0, le=3),
        sorts: str = "Score:desc,MT_IC50_nM:asc,Agretopicity:desc",
        format: str = "tsv",
    ):
        if format not in {"csv", "tsv"}:
            raise HTTPException(422, "Format must be csv or tsv.")
        data = store.read(patient)
        rows = filter_rows(data, q, tier, sorts)
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(
            buffer,
            fieldnames=data["fields"],
            delimiter="\t" if format == "tsv" else ",",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    k: "'" + v
                    if isinstance(v, str) and v.startswith(("=", "+", "-", "@"))
                    else v
                    for k, v in row.items()
                }
            )
        return Response(
            buffer.getvalue(),
            media_type="text/tab-separated-values" if format == "tsv" else "text/csv",
            headers={
                "Content-Disposition": f'attachment; filename="NeoScorer_filtered.{format}"'
            },
        )

    @app.get("/api/inputs")
    def inputs():
        return jobs.inputs()

    @app.post("/api/jobs", status_code=202)
    def start(request: RunRequest):
        return jobs.start(request)

    @app.get("/api/jobs/current")
    def status():
        return jobs.status()

    return app


app = create_app()


def open_when_ready(port):
    for _ in range(100):
        try:
            with urllib.request.urlopen(
                f"http://127.0.0.1:{port}/api/health", timeout=1
            ) as response:
                if json.load(response).get("application") == "neoscorer-dashboard":
                    webbrowser.open(f"http://localhost:{port}")
                    return
        except (OSError, ValueError):
            time.sleep(0.2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8520)
    parser.add_argument("--open-browser", action="store_true")
    args = parser.parse_args()
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as server_socket:
        server_socket.bind(("127.0.0.1", args.port))
        if args.open_browser:
            threading.Thread(
                target=open_when_ready, args=(args.port,), daemon=True
            ).start()
        server = uvicorn.Server(uvicorn.Config(app, host="127.0.0.1", port=args.port))
        server.run(sockets=[server_socket])
