Will be completed once pipeline is completed. 


# Known Limitations

- **Population Δαox** uses a fixed assumed photon index and an i-band UV proxy, without individual Γ filtering, reddening slopes, or Eddington bias correction, giving a larger scatter (~0.33 dex) than fully-cut literature samples (~0.15 dex).
- **Background region contamination** from field-specific sources inside the fixed-radius annulus was confirmed via a 5XMM crossmatch and fixed by excluding a small region around any catalogued source before extraction (`exclude_radius`).
- **Spectral fitting (Γ) had four real bugs, now fixed** — missing background subtraction, a stale post-error-search model state, unpropagated NH/band/redshift/normalisation in MOS combined fits, and boundary-pinned Γ from invalid (negative net-count) detections — with post-fix Γ now matching published values (Nardini et al. 2019).
- **`PhoIndex` hard-limit boundary bug** set the fit's hard limits half a unit beyond the intended `specfit_gamma_min`/`specfit_gamma_max`, letting poorly-constrained fits converge to physically meaningless Γ with a falsely tight error; now fixed to collapse exactly to the intended bounds.
- **Best-fit epoch selection** for multi-ObsID targets was originally by maximum degrees of freedom (sometimes picking the worst-fitting epoch); now by minimum reduced-cstat with a `best_fit_min_dof` floor.
- **Hardness ratio `ZeroDivisionError`** from an unguarded zero source `BACKSCAL` is now caught, with the affected ObsID skipped and logged.
- **Photon index can still be poorly constrained for very low-count sources** — a genuine Poisson-limited floor rather than a code issue, now flagged pre-fit by a configurable `min_counts` gate.
- **MOS-specific Γ–HR decorrelation** (R²≈0.67–0.70 for PN/combined vs. R²≈0.003 for MOS) is consistent with the optical-loading effect reported in Leighly et al. (2007) and isn't correctable in the pipeline; treat MOS-only Γ as less reliable for HR-based conclusions.
- **Fractional variability (Fvar) is unreliable below roughly dof≈200–300**, since the excess-variance estimator's per-bin error term is itself poorly determined at low counts, making the result noise- rather than source-dominated.
- **Δαox's 0.3838 Tananbaum scaling factor was missing** from the internal residual calculation, overstating the metric by ≈2.6×; now fixed.
- **`literature_consistency_check` had four bugs, now fixed** — an error formula that summed absolute confidence bounds instead of taking their width (understating σ by up to ~60×), epoch/fit-type selection that ignored fit reliability and reused an already-fixed max-dof bug pattern, `STACK_GAMMA` fallback errors never being reassigned from the failed `SPEC_GAMMA_PL` attempt, and an unguarded `None` case that could crash the run.
- **Δαox correlates with fit dof** (Spearman ρ≈0.35–0.45, strengthening rather than vanishing after quality filtering), indicating the weak/strong classification isn't fully independent of data quality.
- **An unresolved rest-frame ≈3.5 keV feature** appears episodically in two individually-fit targets, energy-locked but not persistent across epochs, tentatively consistent with K XVIII/Cl XVII rather than anything exotic.
- **Unresolved secondary excess in the population Δαox distribution** at Δαox ≈ −0.6 to −0.8 remains unexplained after ruling out redshift, (g−i) colour, and CIV absorption-index rate as drivers.
- **Cross-catalogue literature comparisons inherit real detection differences** — 4XMM-DR14 and 5XMM-DR15 use independent detection pipelines and stacking, so a non-match doesn't necessarily mean misidentification.
- **Per-target result staleness** is possible since `xray_weakness_comparison` reads whatever per-target files already exist on disk without re-running or version-checking them; confirm results postdate the code version being reported.
- **PN light curves** occasionally fail during `epiclccorr` (a SAS-side segfault on non-standard PI ranges); affected instruments are auto-skipped with `Fvar` left uncomputed.
- **Multi-target parallel runs** could let one run's unconditional `Catalogues/*.fits` deletion race another's read; not a live risk under the current serial-run usage.
