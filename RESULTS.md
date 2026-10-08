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
