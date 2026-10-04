"""Single-worker pipeline jobs with validated inputs and isolated output."""

import json
import os
import subprocess
import sys
import threading
import uuid
from pathlib import Path
from typing import Any, Literal
from fastapi import HTTPException
from pydantic import BaseModel, Field
from patient_io import report_name, validate_inputs


class RunRequest(BaseModel):
    patient_id: str | None = Field(None, max_length=80)
    mutations_file: str | None = None
    expression_file: str | None = None
    mutations: list[dict[str, Any]] | None = None
    expression: dict[str, Any] | None = None
    predictor: Literal["ensemble", "mhcflurry", "netmhcpan", "pwm"] = "ensemble"
    hla: str = Field("HLA-A*02:01", pattern=r"^HLA-[ABC]\*\d{2,3}:\d{2,3}$")
    tpm_threshold: float = Field(1, ge=0, le=1000000, allow_inf_nan=False)
    top: int = Field(12, ge=1, le=1000)


class JobManager:
    def __init__(self, base: Path, reports: Path):
        self.base, self.reports = base, reports
        self.lock = threading.RLock()
        self.job = self.process = self.worker = None
        self.stopping = False

    def inputs(self):
        result = {"mutations": [], "expression": []}
        for path in sorted((self.base / "data").glob("*.json")):
            if path.is_symlink() or path.stat().st_size > 20_000_000:
                continue
            kind = (
                "mutations"
                if "mutation" in path.name.lower()
                else "expression"
                if "expression" in path.name.lower()
                else None
            )
            if kind:
                result[kind].append(path.name)
        return result

    def load_input(self, filename, kind):
        if filename not in self.inputs()[kind]:
            raise ValueError(f"Select an available {kind} input or upload JSON.")
        return json.loads(
            (self.base / "data" / filename).read_text(encoding="utf-8-sig")
        )

    def start(self, request: RunRequest):
        with self.lock:
            if self.stopping or (self.job and self.job["status"] == "running"):
                raise HTTPException(
                    409, "A pipeline job is already running or the server is stopping."
                )
            try:
                mutations = (
                    request.mutations
                    if request.mutations is not None
                    else self.load_input(request.mutations_file, "mutations")
                )
                expression = (
                    request.expression
                    if request.expression is not None
                    else self.load_input(request.expression_file, "expression")
                )
                patient = validate_inputs(
                    mutations, expression, request.patient_id or None
                )
                if request.predictor == "pwm" and request.hla != "HLA-A*02:01":
                    raise ValueError("PWM supports only HLA-A*02:01.")
            except (ValueError, OSError) as exc:
                raise HTTPException(422, str(exc)) from exc
            destination = self.reports / report_name(patient)
            if destination.exists():
                raise HTTPException(
                    409,
                    "This patient already has a report. Use the CLI for an intentional rerun; existing dashboard reports are preserved.",
                )
            job_id = uuid.uuid4().hex
            directory = self.base / ".dashboard" / "jobs" / job_id
            directory.mkdir(parents=True)
            (directory / "mutations.json").write_text(
                json.dumps(mutations), encoding="utf-8"
            )
            (directory / "expression.json").write_text(
                json.dumps(expression), encoding="utf-8"
            )
            self.job = {
                "id": job_id,
                "patient": patient,
                "status": "running",
                "report_id": destination.name.removesuffix(
                    "_clinical_vaccine_prescription.tsv"
                ),
                "message": "Pipeline is starting.",
            }
            command = [
                sys.executable,
                "-u",
                str(self.base / "run_pipeline.py"),
                "--mutations",
                str(directory / "mutations.json"),
                "--expression",
                str(directory / "expression.json"),
                "--patient-id",
                patient,
                "--output-dir",
                str(directory / "output"),
                "--hla",
                request.hla,
                "--predictor",
                request.predictor,
                "--tpm-threshold",
                str(request.tpm_threshold),
                "--top",
                str(request.top),
            ]
            self.worker = threading.Thread(
                target=self.execute, args=(command, directory, destination), daemon=True
            )
            self.worker.start()
            return dict(self.job)

    def execute(self, command, directory, destination):
        try:
            with (directory / "run.log").open("w", encoding="utf-8") as log:
                with self.lock:
                    if self.stopping:
                        raise RuntimeError("Server is shutting down.")
                    self.process = subprocess.Popen(
                        command,
                        cwd=self.base,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        creationflags=subprocess.CREATE_NO_WINDOW
                        if os.name == "nt"
                        else 0,
                    )
                    process = self.process
                code = process.wait(timeout=21600)
            if code != 0:
                raise RuntimeError(
                    f"Pipeline exited with code {code}. See the run log."
                )
            from dashboard_backend.reports import ReportStore

            staged = directory / "output" / destination.name
            ReportStore(staged.parent).read(self.job["report_id"])
            with self.lock:
                if self.stopping:
                    raise RuntimeError("Server stopped before publication.")
                self.reports.mkdir(exist_ok=True)
                # Hard linking publishes a complete file atomically without overwriting.
                os.link(staged.with_suffix(".json"), destination.with_suffix(".json"))
                try:
                    os.link(staged, destination)
                except OSError:
                    destination.with_suffix(".json").unlink(missing_ok=True)
                    raise
                self.job.update(
                    status="completed",
                    message="Report is ready. Patient list refreshed.",
                )
        except Exception as exc:
            with self.lock:
                if self.process and self.process.poll() is None:
                    self.process.kill()
                    self.process.wait()
                if self.job:
                    self.job.update(status="failed", message=str(exc))
        finally:
            with self.lock:
                self.process = None

    def status(self):
        with self.lock:
            if not self.job:
                return {"status": "idle"}
            result = dict(self.job)
            path = self.base / ".dashboard" / "jobs" / result["id"] / "run.log"
            if path.exists():
                with path.open("rb") as stream:
                    stream.seek(max(0, path.stat().st_size - 16000))
                    result["log"] = stream.read().decode("utf-8", errors="replace")
            return result

    def close(self):
        with self.lock:
            self.stopping = True
            if self.process and self.process.poll() is None:
                self.process.terminate()
        if self.worker:
            self.worker.join(timeout=10)
        with self.lock:
            if self.process and self.process.poll() is None:
                self.process.kill()
                self.process.wait()
