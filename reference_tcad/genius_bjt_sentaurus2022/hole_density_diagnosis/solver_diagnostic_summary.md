# Genius NPN BJT hole-density solver diagnosis

## 1. Spatial decomposition

The baseline contains 200 nodes with delta(hQF-psi) between -101 and -99 mV. Of the 204 nodes above 1 decade hole error, 200 are in this plateau.
For SDevice p<1e2 cm^-3, the carrier-relation residual RMSE is 0.00029435 decade.

## 2. Hole continuity terms

Across the plateau, the local hole-row ratio has median 0.150465, P95 1.48757, and maximum 1.97807. The raw Newton hQF correction has median 1.21041 V and is limited to 0.1 V.

## 3. QF update-limit A/B

| Limit | Plateau trial RMSE [decade] | Trial combined residual |
|---:|---:|---:|
| 0.025 V | 1.25902 | 4.26315e-10 |
| 0.05 V | 0.839035 | 4.2631e-10 |
| 0.1 V | 0.00135255 | 4.26312e-10 |
| disabled | 18.1642 | 154769 |

## 4. Strict carrier-row convergence

The enforce/1e-4 run ended after 80 iterations with `carrier_row_convergence`. The original -100 mV plateau was removed and p<1e2 cm^-3 RMSE fell to 0.00228538 decade, but 148 rows still violated the global criterion; the maximum ratio was 0.882285 at node 2672.

## 5. Contact-to-bulk profiles

The base and collector contact values agree. The approximately 0.1 V offset appears only inside the remote n-type low-hole basin, which identifies weak interior minority-carrier coupling rather than a contact boundary-condition mismatch.

## Conclusion

The immediate discrepancy is a prematurely accepted minority-hole QF state. A 0.1 V limited correction removes the plateau and reduces the nonlinear residual; disabling the limit is unstable. Global strict row enforcement is too sensitive to other negligible-density cancellation rows, so the production fix should qualify carrier rows by physical relevance rather than remove the QF update limit.
