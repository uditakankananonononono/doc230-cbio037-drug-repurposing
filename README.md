# doc230-cbio037-drug-repurposing
Audited: no compound-popularity artifact (compound gene-degree AUROC 0.53 vs MLP 0.77); D10 original B1 baseline was degenerate (see RESULTS.md Audit). Verdicts unchanged.
Preregistered benchmark of drug-repurposing link prediction on Hetionet (CPU). See PREREG.md, RESULTS.md, src/run.py. Data: /tmp/het files from hetio/hetionet (hashes in data_manifest). Gates D2 and D3 missed.

## Reproduce / readability audit (2026-10-09)
Data: sources and SHA256 are in data_manifest/. Scripts read fixed paths under /tmp (download the files there first, check the SHA256, then run python3 src/run.py and python3 src/audit.py). Not independently re-run from a clean machine as part of this readability check; the numbers in RESULTS.md come from the original runs on 2 CPU / 2 GB. Hetionet files go in /tmp/het/.
