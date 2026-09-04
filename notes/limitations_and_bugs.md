## Known Limitations.
- **Population Δαox method bias** (fixed Γ assumption, i-band UV proxy) - a deliberate methodological simplification, not a defect.
- **Low-count Γ constraints** - genuine Poisson-statistical floor.
- **MOS-specific Γ–HR decorrelation** - real instrumental physics (optical loading), not something your code got wrong.
- **Fvar reliability threshold below dof≈200–300** - inherent to the excess-variance estimator at low counts.
- **Δαox correlates with fit dof** - a real, important caveat about the weak/strong classification, not a coding error.
- **Unresolved rest-frame ≈3.5 keV feature** - open question, tentatively non-exotic.
- **Unresolved Δαox secondary excess** at -0.6 to -0.8 - open question.
- **Cross-catalogue detection differences** (4XMM-DR14 vs 5XMM-DR15) - a caveat about interpreting literature comparisons, not a defect in your pipeline.
- **Per-target result staleness** - the pipeline behaves exactly as designed (one target per run), the risk is procedural (forgetting to reprocess after a fix), not a bug to fix in code.

## Ongoing Bugs.
- **Multi-target parallel run race condition** - `_Reprocessing()`'s unconditional `Catalogues/*.fits` deletion isn't scoped to the triggering target; a real defect, just not currently triggered because targets are run serially. This is a live bug waiting for a usage pattern that hits it, not something actually fixed.
- **PN light curve `epiclccorr` segfault** - SAS-side, not your code, so not something you can directly fix; currently mitigated by auto-skipping the affected instrument rather than resolved.

## Fixed Bugs.
- **Background region contamination** from field-specific sources inside the annulus - fixed via `exclude_radius`.
- **Spectral fitting (Γ), four bugs** - missing background subtraction, stale post-error-search model state, unpropagated NH/band/redshift/normalisation in MOS combined fits, boundary-pinned Γ from invalid detections.
- **`PhoIndex` hard-limit boundary bug** - hard limits were half a unit beyond the intended config values.
- **Best-fit epoch selection** - was max-dof, now min reduced-cstat with a `best_fit_min_dof` floor.
- **Hardness ratio `ZeroDivisionError`** from unguarded zero source `BACKSCAL`.
- **Δαox's missing 0.3838 Tananbaum scaling factor** - was overstating the metric by ≈2.6×.
- **`literature_consistency_check`, four bugs** - σ formula, epoch/fit-type selection reusing the already-fixed max-dof pattern, `STACK_GAMMA` fallback errors never reassigned, unguarded `None` crash.
