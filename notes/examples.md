# Full Target Results

### Example results using all targets in target_catalogue.json and the parameters set out in config.json:

Γ column shows the best-available fit, preferring combined > PN > MOS (whichever is clean, i.e. not boundary-pinned); "fit" notes which one is shown. HR is the BEHR median with its 68% interval. Fvar: "n.d." = not detected (Fvar² ≤ 0, recorded as upper limit); "—" = not computed.


Found the discrepancy — 6 of your 31 Decent/Good targets aren't in the table you pasted at all (it looks like an older/partial export): `SDSS_J020543.75-063807.0` (Decent strong), `SDSS_J080040.11+110314.1` (Decent weak), `SDSS_J093217.17+285844.8` (Good weak), `SDSS_J122708.29+012638.4` (Decent weak), `SDSS_J131109.57+391128.5` (Good weak), `SDSS_J132040.61+341646.4` (Decent weak). I've left them out below since I have no Γ_own/HR/Fvar for them from your table — let me know if you want me to pull those separately.

For the rest: `Γ_lit` uses `SPEC_GAMMA_PL` (per-obsid power-law fit) when available, falling back to `STACK_GAMMA` (multi-epoch stacked fit) otherwise — I checked both against the 5XMM data and confirmed `SPEC_GAMMA_ERR_LO/UP_PL` are absolute bounds, while `STACK_GAMMA_ERR_LO/UP` are ± deltas (verified `STACK_GAMMA_ERR == mean(LO,UP)` across several rows), so I converted the latter to absolute bounds before use. σ = |Γ_own − Γ_lit| / √(σ_own² + σ_lit²), with each σ taken as half the quoted range — a standard quadrature calc, not a reproduction of your pipeline's exact formula since I still don't have that code, so treat it as a good approximation rather than an exact match.

## Test
| Target | Δαox | Γ (best fit) | Γ_lit | σ | Fit | HR | Fvar |
|---|---|---|---|---|---|---|---|
| 3C 273 | — | 1.70 (1.67, 1.73) | 1.74 (n/a) | — | combined | -0.40 (-0.40, -0.39) | 0.013 ± 0.000 |
| J082619.70+314847.9 | — | 1.79 (1.59, 2.01) | 1.93 (1.66, 2.30) | 0.36 | combined | -0.69 (-1.00, -0.56) | n.d. |
| J094734.19+142116.9 | — | 1.98 (1.92, 2.04) | 1.77 (1.63, 1.95) | 1.20 | combined | -0.69 (-0.74, -0.65) | n.d. |
| J102714.77+354317.4 | — | 1.94 (1.89, 2.00) | 2.10 (1.90, 2.38) | 0.68 | combined | -0.62 (-0.65, -0.59) | 0.153 ± 0.018 |
| PG 1211+143 | — | 2.83 (2.82, 2.83) | 2.36 (2.36, 2.36) | 89.97 | combined | -0.74 (-0.74, -0.74) | 0.130 ± 0.001 |
| PHL 1811 | — | 2.34 (2.26, 2.43) | 2.23 (2.15, 2.30) | 1.01 | combined | -0.55 (-0.60, -0.52) | n.d. |

## Weak (Δαox < -0.3)
| Target | Δαox | Γ (best fit) | Γ_lit | σ | Fit | HR | Fvar |
|---|---|---|---|---|---|---|---|
| J020535.44-043550.8 | -0.38 | 1.82 (1.49, 2.26) | — | — | combined | -0.58 (-0.88, -0.33) | n.d. |
| J021649.79-061152.0 | -0.32 | 1.98 (1.77, 2.21) | 2.01 (1.61, 2.62) | 0.05 | combined | -0.52 (-0.71, -0.34) | n.d. |
| J080342.03+302254.7 | -0.33 | 1.21 (1.14, 1.28) | 1.14 (1.07, 1.24) | 0.60 | combined | -0.21 (-0.27, -0.13) | n.d. |
| J093857.01+412821.2 | -0.35 | 1.09 (0.77, 1.42) | 1.88 (1.62, 2.22) | 1.80 | pn | -0.29 (-0.62, -0.00) | n.d. |
| J114900.13+592225.0 | -0.38 | 0.70 (0.00, 1.48) | 1.31 (1.09, 1.80) | 0.75 | combined | -0.12 (-0.42, 0.22) | — |

## Normal
| Target | Δαox | Γ (best fit) | Γ_lit | σ | Fit | HR | Fvar |
|---|---|---|---|---|---|---|---|
| J022742.83+004002.9 | 0.00 | 1.92 (1.74, 2.11) | 1.85 (1.30, 2.40) | 0.11 | combined | -0.36 (-0.48, -0.24) | — |
| J094446.36+040659.7 | -0.00 | 2.18 (2.01, 2.36) | 2.25 (2.04, 2.51) | 0.22 | combined | -0.82 (-1.00, -0.76) | n.d. |
| J131952.38-022543.3 | -0.11 | 1.73 (1.59, 1.88) | 1.93 (1.70, 2.35) | 0.58 | combined | -0.59 (-0.70, -0.46) | n.d. |
| J142402.13+382104.1 | -0.24 | 1.34 (1.06, 1.63) | 2.75 (2.47, 2.92) | 3.88 | combined | -0.91 (-1.00, -0.82) | n.d. |
| J222406.03-013111.4 | -0.00 | 1.46 (1.23, 1.71) | 2.07 (1.48, 2.67) | 0.95 | combined | -0.24 (-0.53, 0.01) | n.d. |
| J232344.20-005106.9 | -0.01 | 2.10 (1.92, 2.30) | 2.37 (2.21, 2.61) | 0.99 | combined | -0.87 (-1.00, -0.80) | n.d. |
| J234100.91-085531.9 | -0.03 | 2.26 (2.09, 2.44) | 2.08 (1.94, 2.26) | 0.77 | combined | -0.95 (-1.00, -0.91) | n.d. |

## Strong (Δαox > +0.3)
| Target | Δαox | Γ (best fit) | Γ_lit | σ | Fit | HR | Fvar |
|---|---|---|---|---|---|---|---|
| J002239.10+012950.2 | 0.34 | 2.07 (1.99, 2.15) | 1.82 (1.68, 2.00) | 1.42 | combined | -0.66 (-0.71, -0.61) | n.d. |
| J085141.76+161221.9 | 0.32 | 1.23 (0.89, 1.63) | 1.93 (1.76, 2.15) | 1.67 | combined | -0.11 (-0.28, 0.08) | — |
| J090404.15+151254.5 | 0.37 | 1.94 (1.90, 1.98) | 1.92 (1.90, 1.95) | 0.33 | combined | -0.61 (-0.64, -0.58) | n.d. |
| J092541.34+131938.9 | 0.32 | 2.17 (2.09, 2.25) | 2.07 (1.94, 2.23) | 0.60 | combined | -0.69 (-0.74, -0.64) | n.d. |
| J123508.21+391419.9 | 0.36 | 1.88 (1.72, 2.06) | 1.82 (1.63, 2.02) | 0.24 | combined | -0.62 (-0.75, -0.47) | n.d. |
| J125005.72+263107.5 | 0.60 | 2.08 (2.05, 2.10) | 2.08 (2.05, 2.11) | 0.04 | combined | -0.63 (-0.64, -0.61) | n.d. |
| J212951.16+004808.8 | 0.36 | 1.87 (1.81, 1.92) | 1.84 (1.74, 1.97) | 0.25 | combined | -0.57 (-0.61, -0.53) | n.d. |
# Population Summary

### Example results of 4XMM-DR14 and SDSS DR16Q using the parameters set out in config.json:

4XMM-DR14 × SDSS DR16Q crossmatch, z = 1.8–2.2, after quality cuts (SC_SUM_FLAG, PSFFLUX, colour, radio, BAL).

## Sample composition

| | N | Fraction |
|---|---|---|
| Total | 1539 | — |
| Weak (Δαox < -0.3) | 269 | 17.5% |
| Normal (-0.3 ≤ Δαox ≤ 0.3) | 1005 | 65.3% |
| Strong (Δαox > +0.3) | 265 | 17.2% |

## Δαox distribution

| Statistic | Value |
|---|---|
| Standard deviation | 0.333 dex |
| Mean | ≈ 0.00 (consistent with zero) |

## αox

| Statistic | Value |
|---|---|
| Median | -1.641 |
| 16th percentile | -1.802 |
| 84th percentile | -1.499 |

## Luminosities (population median)

| Quantity | Value |
|---|---|
| log L(2 keV) | 26.39 erg s⁻¹ Hz⁻¹ |
| log L(2500 Å) | 30.63 erg s⁻¹ Hz⁻¹ |

## Classification thresholds

| Parameter | Value |
|---|---|
| xray_weak_threshold | -0.3 |
| xray_strong_threshold | +0.3 |
| z_min – z_max | 1.8 – 2.2 |
