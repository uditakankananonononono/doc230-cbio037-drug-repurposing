# CBIO037 results (results/results.json from src/run.py; log results/run.log). 5 seeds, mean AUROC (sd in json). Hetionet v1.0, CtD = 755 edges, 387 compounds, 77 diseases.
| Setting | B1 degree | B2 LR | B3 gene cosine | M MLP |
|---|---|---|---|---|
| cold-compound (D2) | 0.855 | 0.738 | 0.586 | 0.766 |
| random split (D3) | 0.915 | 0.890 | 0.542 | 0.893 |
| cold-disease (D4) | 0.738 | 0.711 | 0.574 | 0.732 |
| label-shuffle control (D5) | - | 0.506 | - | 0.506 |
| degree-matched negatives, cold-compound (D6) | 0.587 | 0.616 | 0.578 | 0.725 |
| degree-matched negatives, random (D6) | 0.655 | 0.722 | 0.567 | 0.810 |
| gene profiles only (D7) | - | 0.755 | - | 0.775 |
| side-effect + disease profile (D7, see caveat) | - | 0.762 | - | 0.800 |
| CpD target, cold-compound (D8) | 0.894 | 0.728 | 0.669 | 0.740 |
D9: seed sd 0.014-0.07, listed in results.json. D10 hits@10 per held-out disease: in results.json.
Gates: D2 NOT met (M 0.766 vs degree baseline 0.855 with random negatives). D3 NOT met (random 0.893 minus cold 0.766 = 0.127 > 0.05): performance drops on unseen compounds.
Honest reading: with random negatives, a degree-product baseline beats the MLP; most apparent signal is node popularity. With degree-matched negatives the degree baseline falls to 0.59 and M keeps 0.725 (cold-compound), so the MLP has signal beyond degree, but it is modest. Shuffle control is at chance (0.506), so no pipeline leakage found.
Caveats: D7 "side-effect only" also includes the disease gene profile, so it is not a pure side-effect test. Hetionet CtD edges were curated from older sources; negatives are unlabelled pairs, not known non-treatments. Top-ranked predictions are not validated repurposing candidates. No GNN or metapath methods run. Replication/benchmark, no novelty claim.
D10 hits@10 (fraction of held-out disease instances, n=252 over 5 seeds, candidates = held-out CtD compounds + 200 non-CtD): M 0.377, B1 0.115, B3 0.151. M is best at ranking here.

## Audit (second pass) - no compound-popularity artifact found; D10 baseline was degenerate; verdicts unchanged
Prereg AUDIT_PREREG.md committed before running; src/audit.py, results/audit.json. Data SHA256 re-verified; original D2 M AUROC reproduced exactly (0.766).
Correction to my own prereg premise: in the cold-compound split the TEST negatives were already drawn only from the held-out positive compounds (the code sets the test compound pool to the held-out set with CtD), so the test set was already compound-matched. That is why B3 and B4 AUROC are identical in both audit modes below (same test pairs); only the TRAINING negatives differ between modes.
- Q1 B4 (disease-independent compound gene-degree) cold-compound AUROC 0.528 (sd 0.047) vs M 0.766: B4 << M. The pre-specified "annotation density explains M" condition is NOT met. No compound-popularity artifact in the test AUROC.
- Q2 Training negatives drawn from the positives' own compound pool (compound-matched): M 0.832 (sd 0.024), B2 0.793, B3 0.586, B4 0.528. M > 0.60 bar met: a disease-specific signal exists on cold compounds (and it is stronger when training negatives are compound-matched). Note original D2 gate stays NOT met because the disease-degree baseline B1 (0.855 with random negatives) beats M; that comparison is about disease popularity, which Q1/Q2 do not touch (D6 degree-matched negatives remain the control: B1 0.587, M 0.725).
- Q3 D10 hits@10 (252 held-out disease instances, 5 seeds): M 0.377, B3 0.151, B4 (compound gene-degree) 0.270, random expectation 0.089. The original B1 baseline (0.115) was degenerate (all candidates have train CtD degree 0, ranked by random tie-break), so "M best at ranking" was against a weak baseline. Against the disease-independent B4, M still leads by 0.107 (> 0.05 bar met), but B4's 0.27 vs random 0.089 shows the candidate set itself leaks popularity (held-out CtD compounds are better annotated than the random non-CtD ones), so D10 hits@10 is inflated for every method.
Revised statement: M (MLP on gene-profile features) has a real but modest disease-specific signal on unseen compounds (0.73 with degree-matched negatives, 0.83 with compound-matched training negatives); D2 NOT met and D3 NOT met stand. D10 should be cited only with the B4 comparison. Not validated repurposing candidates.
