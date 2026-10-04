# NeoScorer multi-patient dashboard

A local FastAPI application with a vanilla JavaScript frontend and bundled Chart.js 4.4.8. Patient discovery, report parsing, job execution and visualization are separate modules. No frontend build or internet connection is required to explore existing reports.

## Start

Double-click `C:\Tools\LAUNCHERS\Launch_NeoScorer_Dashboard.bat`, or run:

```powershell
& "C:\Tools\.venv\Scripts\python.exe" "C:\Tools\Mini-NeoScorer\dashboard.py" --open-browser
```

The app binds only to `127.0.0.1:8520`. Use `--port 8522` for a different port. The browser opens after the health endpoint responds. Ctrl+C stops the server and its active pipeline worker. A port conflict fails without terminating other applications. Run one dashboard instance against a reports directory; there is one active job slot per server.

## Patient reports and evidence

The backend rescans `reports/*_clinical_vaccine_prescription.tsv` on every discovery request. Add a report and click the refresh button, or complete a new run. No patient ID or filename is embedded in the dashboard. The filename prefix is the discovery key; underscores display as hyphens. New pipeline outputs additionally preserve the original patient ID in a JSON sidecar.

- All original prescription columns are shown. Model-specific fields appear in the selected-candidate card and remain in exports.
- Alterations screened means **rows in the report**, not the original WES count. The pipeline exports one best epitope per successfully screened variant.
- Strong binders satisfy **positive IC50 ≤ 50 nM OR percentile between 0 and 0.5%**. The smaller affinity-only count is also shown.
- Tier counts use recorded `ClinicalTier` assignments, never an invented payload size.
- HLA comes from `HLA`, `HLA_Allele` or `hla_allele`, then sidecar metadata. Legacy reports fall back to HLA references in `Action`, explicitly labeled as inferred. Missing HLA is shown as not recorded.
- Source WES counts, candidate 9-mer totals and RNA gates are displayed only when recorded in the sidecar. No historical or patient-specific constants are substituted.
- The quadrant guides are exploratory reference thresholds: TPM 5 and IC50 50 nM. Non-positive or missing values are omitted from logarithmic charts with an omission count. Unknown/missing metrics remain missing.
- MHCflurry presentation is an integrated score. The ensemble rank is a combined rank and is not presented as a separate NetMHCpan rank. `MHCflurry_Filter_NonBinder` sentinel values are not shown as genuine NetMHCpan predictions. Individual model ranks and isolated cleavage probabilities are not stored by the existing pipeline.
- These are research predictions, with no added clinical validation or treatment recommendations.

## Exploration and exports

Search gene, mutation, mutant peptide or wild-type peptide. Tier filtering updates both charts and the grid; the top-level patient KPIs remain whole-report totals. Click a table heading to sort, or Shift-click to add a sorting level. The sorting dropdown supplies score, affinity and agretopicity ordering with deterministic secondary keys. Click a gene or chart point/bar to inspect model evidence. Tables paginate at 20 rows while charts and exports use all matching rows.

CSV and TSV downloads preserve original columns and the active search, tier filter and complete sorting order. Spreadsheet formula prefixes are escaped in text fields. Patient switching resets filters and selection; stale asynchronous responses cannot overwrite a newer patient's view.

## New patients

Open **Run new patient** and select JSON files discovered in `data/`, or upload files from your computer. Uploaded JSON takes precedence over a local selection; files are validated before the worker starts. The request limit is 20 MB. A run stores its inputs, console log and staged output under the git-ignored `.dashboard/jobs/<job-id>/` directory.

Mutation input is a non-empty array of cBioPortal-compatible objects. Each entry needs `gene.hugoGeneSymbol`, `entrezGeneId`, `mutationType`, `tumorAltCount`, and `tumorRefCount`; missense entries also need `proteinChange`. Expression input is a non-empty object mapping Entrez gene IDs to finite, non-negative TPM values. The patient ID is inferred from mutation `patientId` or entered explicitly. Mixed-patient inputs and identity mismatches are rejected. Expression JSON has no patient identity field, so the user must supply the matched RNA sample.

The sidebar controls predictor, HLA, TPM threshold and the number of top candidates printed in the run log. **Top** does not truncate the saved report. Each run is a subprocess of the shared Python environment. Status/logs update while the page is open; completion discovers and selects the new patient. Failed jobs retain their logs without publishing a report. Jobs time out after six hours. Status is in memory; input snapshots and logs persist across server restarts, but interrupted jobs do not resume automatically.

Existing patient reports are never overwritten by the web interface. For an intentional rerun, use the CLI:

```powershell
& "C:\Tools\.venv\Scripts\python.exe" "C:\Tools\Mini-NeoScorer\run_pipeline.py" `
  --mutations "C:\path\patient_mutations.json" `
  --expression "C:\path\patient_expression.json" `
  --patient-id "PATIENT-123" --hla "HLA-A*02:01" `
  --predictor ensemble --tpm-threshold 1 --top 12
```

The CLI now supports `--mutations`, `--expression`, `--patient-id`, and `--output-dir` in addition to its existing arguments. Existing default inputs still work, but output names follow the input patient identity. New reports include `HLA_Allele` and a matching JSON sidecar with patient identity, HLA, predictor, threshold, timestamp and funnel statistics.

MHCflurry requires its installed model weights. Ensemble and NetMHCpan use NIH IEDB and transmit peptide sequences and HLA. Protein resolution can query UniProt. These existing pipeline behaviors are disclosed in the run panel. PWM supports only HLA-A*02:01 and is labeled a baseline. No patient inputs are sent by the dashboard merely for report viewing. Model inference behavior itself was not redesigned in this change.

## Layout

- `dashboard.py`: application factory, API routes, origin/host guards, local entrypoint.
- `dashboard_backend/reports.py`: report discovery, normalization, HLA evidence, filtering/sorting.
- `dashboard_backend/jobs.py`: validated run requests, single worker, logs and output publication.
- `patient_io.py`: shared identity, schema validation and report naming.
- `dashboard_static/index.html`, `styles.css`: responsive dark interface.
- `dashboard_static/data.js`: API access, filtering, sorting and display helpers.
- `dashboard_static/charts.js`: chart lifecycle and logarithmic visualizations.
- `dashboard_static/app.js`: patient state, grid, consensus and run-panel interactions.
- `dashboard_static/vendor/`: pinned local Chart.js bundle and MIT license.
- `tests/test_dashboard.py`: isolated API, pipeline and patient-isolation regression tests.

This is a local single-operator application, not a network-deployed clinical service with authentication, audit infrastructure or validated clinical decision support. Do not expose the server on a public interface. Uploaded input snapshots and logs persist locally for troubleshooting; manage their retention according to your workspace's data policy.

## Verification

```powershell
& "C:\Tools\.venv\Scripts\python.exe" -m pytest tests\test_dashboard.py -q
& "C:\Tools\.venv\Scripts\python.exe" -m ruff check dashboard.py patient_io.py dashboard_backend tests\test_dashboard.py run_pipeline.py --select E9,F
& "C:\Tools\.venv\Scripts\python.exe" -m compileall -q dashboard.py patient_io.py dashboard_backend run_pipeline.py tests
node --check dashboard_static/app.js
node --check dashboard_static/charts.js
node --check dashboard_static/data.js
```

Regression tests use temporary directories and include a real subprocess run with a synthetic silent variant, avoiding remote prediction calls or changes to the patient's report. Browser checks cover search, empty results, tier selection, sorting, CSV export, consensus selection, patient switching, new-run completion and desktop/mobile layouts. Full ensemble model inference is not part of the dashboard smoke tests.
