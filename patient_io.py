"""Patient identity, JSON validation and report naming shared by CLI and dashboard."""

import math
import re

REPORT_SUFFIX = "_clinical_vaccine_prescription.tsv"


def validate_inputs(mutations, expression, patient_id=None):
    if (
        not isinstance(mutations, list)
        or not mutations
        or not all(isinstance(m, dict) for m in mutations)
    ):
        raise ValueError("Mutations must be a non-empty JSON array of variant objects.")
    if not isinstance(expression, dict) or not expression:
        raise ValueError(
            "Expression must be a non-empty JSON object mapping Entrez gene IDs to TPM."
        )
    if any(
        isinstance(v, bool)
        or not isinstance(v, (int, float))
        or not math.isfinite(v)
        or v < 0
        for v in expression.values()
    ):
        raise ValueError("Expression values must be finite, non-negative numbers.")
    identities = {str(m["patientId"]) for m in mutations if m.get("patientId")}
    if len(identities) > 1:
        raise ValueError(
            "Mutation input contains multiple patients. Split it into individual patient inputs."
        )
    inferred = next(iter(identities), None)
    patient_id = patient_id or inferred
    if not patient_id or not re.fullmatch(
        r"[A-Za-z0-9][A-Za-z0-9_-]{0,79}", patient_id
    ):
        raise ValueError(
            "Provide a patient ID using 1–80 letters, digits, hyphens or underscores."
        )
    if inferred and patient_id != inferred:
        raise ValueError(f"Patient ID does not match the mutation input ({inferred}).")
    for m in mutations:
        if (
            not isinstance(m.get("gene"), dict)
            or not m["gene"].get("hugoGeneSymbol")
            or "entrezGeneId" not in m
            or not isinstance(m.get("mutationType"), str)
        ):
            raise ValueError(
                "Each variant needs gene.hugoGeneSymbol, entrezGeneId and mutationType (cBioPortal format)."
            )
        for key in ("tumorAltCount", "tumorRefCount"):
            v = m.get(key)
            if isinstance(v, bool) or not isinstance(v, int) or v < 0:
                raise ValueError(f"Each variant needs a non-negative integer {key}.")
        if m["mutationType"] == "Missense_Mutation" and not isinstance(
            m.get("proteinChange"), str
        ):
            raise ValueError("Missense variants need proteinChange.")
    return patient_id


def report_name(patient_id):
    return patient_id.replace("-", "_") + REPORT_SUFFIX
