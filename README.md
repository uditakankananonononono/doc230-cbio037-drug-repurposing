negative (audited): preregistered Hetionet drug-repurposing link prediction missed gates D2 and D3; no compound-popularity artifact (gene-degree AUROC 0.53 vs MLP 0.77), but the D10 hits@10 headline (M 0.377) is inflated: the original B1 baseline was degenerate, and per the RESULTS.md Audit, "B4's 0.27 vs random 0.089 shows the candidate set itself leaks popularity (held-out CtD compounds are better annotated than the random non-CtD ones), so D10 hits@10 is inflated for every method" (builder-reported).

# doc230-cbio037-drug-repurposing
Preregistered benchmark of drug-repurposing link prediction on Hetionet (CPU). See PREREG.md, RESULTS.md, src/run.py. Data: /tmp/het files from hetio/hetionet (hashes in data_manifest). Gates D2 and D3 missed.

## Reproduce / readability audit (2026-10-09)
Data: sources and SHA256 are in data_manifest/. Scripts read fixed paths under /tmp (download the files there first, check the SHA256, then run python3 src/run.py and python3 src/audit.py). Not independently re-run from a clean machine as part of this readability check; the numbers in RESULTS.md come from the original runs on 2 CPU / 2 GB. Hetionet files go in /tmp/het/.
