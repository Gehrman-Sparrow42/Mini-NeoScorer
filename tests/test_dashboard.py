"""Regression coverage for patient isolation, exports, inputs and pipeline jobs."""

import csv
import io
import json
import shutil
import time
from unittest.mock import patch
import pytest
from fastapi.testclient import TestClient
from dashboard import BASE, create_app
from patient_io import validate_inputs

SAMPLE = {
    "Rank": "1",
    "Gene": "TARGET_A",
    "ProteinChange": "A2V",
    "Neopeptide_9mer": "AVAAAAAAA",
    "WT_9mer": "AAAAAAAAA",
    "MutPos": "2",
    "AF": "0.5",
    "Reads": "10/20",
    "TPM": "10",
    "MT_IC50_nM": "20",
    "Percentile_Rank": "0.1",
    "Presentation_Score": "0.8",
    "MHCflurry_IC50": "25",
    "NetMHCpan_IC50": "16",
    "WT_IC50_nM": "100",
    "Agretopicity": "5",
    "Score": "90",
    "ClinicalTier": "Tier 1: High Priority (Vaccine Payload)",
    "Action": "INCLUDE: HLA-A*02:01 presentation.",
    "Predictor": "Consensus_Ensemble(MHCflurry+NetMHCpan)",
}
SAMPLE_ROWS = [
    SAMPLE,
    {
        **SAMPLE,
        "Rank": "2",
        "Gene": "BACKUP_B",
        "MT_IC50_nM": "100",
        "Percentile_Rank": "0.4",
        "Score": "80",
        "ClinicalTier": "Tier 2: Secondary / Backup",
        "Neopeptide_9mer": "ABAAAAAAA",
    },
    {
        **SAMPLE,
        "Rank": "3",
        "Gene": "ANCHOR_C",
        "MT_IC50_nM": "40",
        "Percentile_Rank": "0.3",
        "Agretopicity": "20",
        "Neopeptide_9mer": "ACAAAAAAA",
    },
    {
        **SAMPLE,
        "Rank": "4",
        "Gene": "EXCLUDED_D",
        "MT_IC50_nM": "5000",
        "Percentile_Rank": "50",
        "Score": "70",
        "ClinicalTier": "Tier 3: Deprioritized",
        "Neopeptide_9mer": "ADAAAAAAA",
    },
]

HEADERS = {"X-NeoScorer": "dashboard"}


@pytest.fixture
def workspace(tmp_path):
    (tmp_path / "reports").mkdir()
    (tmp_path / "data").mkdir()
    with (tmp_path / "reports" / "PATIENT_A_clinical_vaccine_prescription.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(SAMPLE), delimiter="\t")
        writer.writeheader()
        writer.writerows(SAMPLE_ROWS)
    rows = SAMPLE_ROWS
    with (tmp_path / "reports" / "PATIENT_B_clinical_vaccine_prescription.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(
            stream, fieldnames=[*rows[0], "HLA_Allele"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerow({**rows[0], "Gene": "=FORMULA", "HLA_Allele": "HLA-B*07:02"})
    return tmp_path


@pytest.fixture
def client(workspace):
    with TestClient(create_app(workspace)) as client:
        yield client


def test_patient_isolation_and_hla(client):
    patients = client.get("/api/patients").json()
    assert [p["patient"] for p in patients] == ["PATIENT-A", "PATIENT-B"]
    a = client.get("/api/patients/PATIENT_A").json()
    b = client.get("/api/patients/PATIENT_B").json()
    assert a["total"] == 4 and b["total"] == 1
    assert a["tier1"] == 2 and a["strong_affinity"] == 2 and a["strong"] == 3
    assert a["hla_source"].startswith("Inferred")
    assert b["hla"] == "HLA-B*07:02" and b["hla_source"] == "Report column"
    assert client.get("/api/patients/nonexistent").status_code == 404
    assert client.get("/api/patients/..%2Fdata").status_code == 404


def test_exports_filter_sort_and_escape(client):
    response = client.get(
        "/api/patients/PATIENT_A/export",
        params={"tier": 1, "sorts": "Agretopicity:desc,Score:desc", "format": "csv"},
    )
    rows = list(csv.DictReader(io.StringIO(response.text)))
    assert len(rows) == 2
    assert rows[0]["Gene"] == "ANCHOR_C"
    assert all(r["ClinicalTier"].startswith("Tier 1") for r in rows)
    assert "attachment" in response.headers["content-disposition"]
    result = client.get(
        "/api/patients/PATIENT_A/export", params={"q": "avaaaaaaa"}
    ).text
    assert len(list(csv.DictReader(io.StringIO(result), delimiter="\t"))) == 1
    assert "'=FORMULA" in client.get("/api/patients/PATIENT_B/export").text
    assert client.get("/api/patients/PATIENT_A/export?sorts=bad:asc").status_code == 422
    assert client.get("/api/patients/PATIENT_A/export?format=html").status_code == 422


def test_empty_and_malformed_reports(client, workspace):
    fields = "\t".join(SAMPLE)
    (workspace / "reports" / "EMPTY_clinical_vaccine_prescription.tsv").write_text(
        fields + "\n"
    )
    empty = client.get("/api/patients/EMPTY").json()
    assert empty["total"] == 0 and empty["hla"] == "Not recorded"
    (workspace / "reports" / "BAD_clinical_vaccine_prescription.tsv").write_text(
        "broken\n"
    )
    assert client.get("/api/patients/BAD").status_code == 422
    assert client.get("/api/health").status_code == 200


def test_origin_and_validation(client):
    assert client.post("/api/jobs", json={}).status_code == 403
    assert (
        client.post(
            "/api/jobs",
            headers={**HEADERS, "Origin": "https://untrusted.invalid"},
            json={},
        ).status_code
        == 403
    )
    assert (
        client.post(
            "/api/jobs", headers=HEADERS, json={"mutations_file": "../secret.json"}
        ).status_code
        == 422
    )
    assert (
        client.post(
            "/api/jobs", headers=HEADERS, json={"tpm_threshold": -1}
        ).status_code
        == 422
    )
    assert (
        client.get("/api/health", headers={"host": "untrusted.invalid"}).status_code
        == 400
    )


def inputs(patient="TEST-NEW"):
    return [
        {
            "patientId": patient,
            "gene": {"hugoGeneSymbol": "SYNTHETIC"},
            "entrezGeneId": 1,
            "mutationType": "Silent",
            "tumorAltCount": 1,
            "tumorRefCount": 2,
        }
    ], {"1": 10}


def test_validation_patient_identity():
    mutations, expression = inputs()
    assert validate_inputs(mutations, expression) == "TEST-NEW"
    with pytest.raises(ValueError):
        validate_inputs(mutations, expression, "OTHER")
    with pytest.raises(ValueError):
        validate_inputs(mutations, {"1": float("nan")})
    with pytest.raises(ValueError):
        validate_inputs(mutations + inputs("OTHER")[0], expression)


def test_existing_report_not_overwritten(client):
    mutations, expression = inputs("PATIENT-A")
    response = client.post(
        "/api/jobs",
        headers=HEADERS,
        json={"mutations": mutations, "expression": expression, "predictor": "pwm"},
    )
    assert response.status_code == 409


def test_real_background_cli_job(client, workspace):
    # Silent synthetic variant exercises the real subprocess/CLI/report publication
    # without submitting peptide sequences or invoking model downloads.
    shutil.copy(BASE / "run_pipeline.py", workspace / "run_pipeline.py")
    shutil.copy(BASE / "patient_io.py", workspace / "patient_io.py")
    shutil.copytree(
        BASE / "src", workspace / "src", ignore=shutil.ignore_patterns("__pycache__")
    )
    mutations, expression = inputs()
    response = client.post(
        "/api/jobs",
        headers=HEADERS,
        json={"mutations": mutations, "expression": expression, "predictor": "pwm"},
    )
    assert response.status_code == 202, response.text
    for _ in range(100):
        state = client.get("/api/jobs/current").json()
        if state["status"] != "running":
            break
        time.sleep(0.1)
    assert state["status"] == "completed", state
    report = client.get("/api/patients/TEST_NEW").json()
    assert report["total"] == 0 and report["hla"] == "HLA-A*02:01"
    assert report["metadata"]["stats"]["total_somatic_mutations"] == 1
    assert client.get("/api/patients/PATIENT_A").json()["total"] == 4


def test_failed_background_run(client, workspace):
    (workspace / "run_pipeline.py").write_text(
        'raise RuntimeError("deliberate test failure")\n'
    )
    mutations, expression = inputs()
    response = client.post(
        "/api/jobs",
        headers=HEADERS,
        json={"mutations": mutations, "expression": expression, "predictor": "pwm"},
    )
    assert response.status_code == 202
    for _ in range(50):
        state = client.get("/api/jobs/current").json()
        if state["status"] != "running":
            break
        time.sleep(0.05)
    assert state["status"] == "failed" and "deliberate test failure" in state["log"]
    assert client.get("/api/patients/TEST_NEW").status_code == 404


def test_pipeline_hla_and_identity(tmp_path):
    import run_pipeline

    mutations, expression = inputs("NEW-PATIENT")
    m, e = tmp_path / "m.json", tmp_path / "e.json"
    m.write_text(json.dumps(mutations))
    e.write_text(json.dumps(expression))
    with patch.object(
        run_pipeline,
        "screen_patient_mutations",
        return_value=(
            [],
            dict.fromkeys(
                [
                    "total_somatic_mutations",
                    "missense_mutations",
                    "gated_out_low_rna",
                    "isoform_mismatch_or_unresolved",
                    "successfully_screened_variants",
                    "total_9mers_evaluated",
                    "strong_binders_found",
                    "weak_binders_found",
                    "non_binders_found",
                ],
                0,
            ),
        ),
    ):
        run_pipeline.run_pipeline(
            m,
            e,
            hla_allele="HLA-B*07:02",
            patient_id="NEW-PATIENT",
            output_dir=tmp_path,
        )
    metadata = json.loads(
        (tmp_path / "NEW_PATIENT_clinical_vaccine_prescription.json").read_text()
    )
    assert (
        metadata["patient_id"] == "NEW-PATIENT"
        and metadata["hla_allele"] == "HLA-B*07:02"
    )
    assert (
        "HLA_Allele"
        in (tmp_path / "NEW_PATIENT_clinical_vaccine_prescription.tsv").read_text()
    )
