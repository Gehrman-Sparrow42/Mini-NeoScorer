"""Filesystem-backed reports; no patient-specific constants."""

import csv
import json
import math
import re
from pathlib import Path

from fastapi import HTTPException
from patient_io import REPORT_SUFFIX

NUMERIC = {
    "Rank",
    "MutPos",
    "AF",
    "TPM",
    "MT_IC50_nM",
    "Percentile_Rank",
    "Presentation_Score",
    "MHCflurry_IC50",
    "NetMHCpan_IC50",
    "WT_IC50_nM",
    "Agretopicity",
    "Score",
}
REQUIRED = {
    "Rank",
    "Gene",
    "ProteinChange",
    "Neopeptide_9mer",
    "WT_9mer",
    "MutPos",
    "AF",
    "Reads",
    "TPM",
    "MT_IC50_nM",
    "Percentile_Rank",
    "Agretopicity",
    "Score",
    "ClinicalTier",
    "Action",
    "Predictor",
}


class ReportStore:
    def __init__(self, directory: Path):
        self.directory = directory

    def discover(self):
        return [
            {
                "id": p.name[: -len(REPORT_SUFFIX)],
                "patient": p.name[: -len(REPORT_SUFFIX)].replace("_", "-"),
                "file": p.name,
            }
            for p in sorted(self.directory.glob("*" + REPORT_SUFFIX))
            if p.is_file() and not p.is_symlink()
        ]

    def path(self, key):
        match = next((p for p in self.discover() if p["id"] == key), None)
        if match is None:
            raise HTTPException(404, "Patient report not found.")
        return self.directory / match["file"]

    def read(self, key):
        path = self.path(key)
        try:
            with path.open(encoding="utf-8-sig", newline="") as stream:
                reader = csv.DictReader(stream, delimiter="\t")
                fields = reader.fieldnames or []
                if not REQUIRED.issubset(fields):
                    raise ValueError(
                        "Missing columns: " + ", ".join(sorted(REQUIRED - set(fields)))
                    )
                rows = []
                for index, raw in enumerate(reader):
                    if None in raw or any(v is None for v in raw.values()):
                        raise ValueError(f"Malformed TSV row {index + 2}.")
                    row = dict(raw)
                    for field in NUMERIC:
                        try:
                            val = float(raw.get(field) or "nan")
                            row[field] = val if math.isfinite(val) else None
                        except ValueError:
                            row[field] = None
                    match = re.match(
                        r"Tier\s+([123])\b", raw.get("ClinicalTier", ""), re.I
                    )
                    row["tier"] = int(match[1]) if match else 0
                    row["uid"] = index
                    rows.append(row)
            metadata = {}
            sidecar = path.with_suffix(".json")
            if sidecar.exists():
                metadata = json.loads(sidecar.read_text(encoding="utf-8"))
                if not isinstance(metadata, dict):
                    raise ValueError("Metadata must be a JSON object.")
            hla = {
                str(r[k])
                for r in rows
                for k in ("HLA", "HLA_Allele", "hla_allele")
                if r.get(k)
            }
            source = "Report column"
            if not hla and metadata.get("hla_allele"):
                hla = {str(metadata["hla_allele"])}
                source = "Run metadata"
            if not hla:
                hla = {
                    m
                    for r in rows
                    for m in re.findall(r"HLA-[A-Z]+\*\d+:\d+", r.get("Action", ""))
                }
                source = "Inferred from report Action text" if hla else "Not recorded"
            strong_affinity = sum(
                r["MT_IC50_nM"] is not None and 0 < r["MT_IC50_nM"] <= 50 for r in rows
            )
            strong = sum(
                (r["MT_IC50_nM"] is not None and 0 < r["MT_IC50_nM"] <= 50)
                or (
                    r["Percentile_Rank"] is not None
                    and 0 <= r["Percentile_Rank"] <= 0.5
                )
                for r in rows
            )
            return {
                "id": key,
                "patient": metadata.get("patient_id") or key.replace("_", "-"),
                "source": path.name,
                "modified": path.stat().st_mtime,
                "hla": ", ".join(sorted(hla)) or "Not recorded",
                "hla_source": source,
                "rows": rows,
                "fields": fields,
                "total": len(rows),
                "strong": strong,
                "strong_affinity": strong_affinity,
                "tier1": sum(r["tier"] == 1 for r in rows),
                "metadata": metadata,
            }
        except (OSError, ValueError, csv.Error) as exc:
            raise HTTPException(422, f"Cannot read report: {exc}") from exc


def filter_rows(
    report, query="", tier=0, sorts="Score:desc,MT_IC50_nM:asc,Agretopicity:desc"
):
    rows = [
        r
        for r in report["rows"]
        if (not tier or r["tier"] == tier)
        and query.casefold()
        in " ".join(
            str(r.get(k, ""))
            for k in ("Gene", "ProteinChange", "Neopeptide_9mer", "WT_9mer")
        ).casefold()
    ]
    specs = []
    for item in sorts.split(","):
        parts = item.split(":")
        if (
            len(parts) != 2
            or parts[0] not in REQUIRED | NUMERIC
            or parts[1] not in {"asc", "desc"}
        ):
            raise HTTPException(422, "Invalid sort specification.")
        specs.append(parts)
    for key, direction in reversed(specs):
        rows = sorted(
            [r for r in rows if r.get(key) is not None],
            key=lambda r: r[key],
            reverse=direction == "desc",
        ) + [r for r in rows if r.get(key) is None]
    return rows
