# X-ray Strength AGN Code
### NOTE: Project still under development, README may not be up-to-date. Code subject to change and may contain bugs.

Modular pipeline created to analyse X-ray weak, normal and strong AGN from XMM-Newton observations.
HEASoft and SAS must be installed and initialised (sourced) before running. BEHR must be compiled and its binary placed in the project root directory (not on PATH).

Initiate with `chmod +x run_pipeline.sh`.
Run program as `./run_pipeline.sh` in BASH.

---

## Config parameters

### Target & working directory
- **target_name**: Name of the target you wish to analyse.
- **working_dir**: Directory for the program to operate in. Use `"working"` for the directory containing the `.py` file, or provide an absolute path.
- **field**: Cosmology you wish to use in analysis. Options can be found in `field_catalogue.json`.

### Download
- **ObsID_to_exclude**: ObsIDs to skip. Example: `["0060010101", "0060010102"]`.
- **ObsID_specifically_chosen**: ObsIDs to use (downloads only these). Example: `["0060010101", "0060010102"]`.
- **number_of_datasets**: How many datasets wanted in total. If fewer exist, the program downloads what's available.
- **concurrent_downloads**: Number of datasets to download concurrently. Recommended 5-10.
- **CCF_path**: Path to CCF folder.

### General
- **reprocess**: `false` to skip reprocessing entirely. `true` to regenerate files even if they already exist, with no confirmation. Any other string (e.g. `"careful"`) also regenerates, but first prompts `[y/N]` in the terminal before deleting anything, since this step permanently removes existing processed files — useful when you want the safety net but not the silence of `false`, or the risk of `true`.
- **auto_resolve**: `true` to backtrack and run missing prerequisites automatically.
- **min_counts**: Minimum counts to prevent running whole pipeline for target too faint to obtain a photon index value.
- **stop_at_low_counts**:  `true` to stop pipeline if number of counts of a target is lower than the `min_counts` value.

### Population crossmatch
- **crossmatch_radius_arcsec**: Matching radius used when crossmatching 4XMM against SDSS DR16.
- **z_min**: Lower redshift bound for the population sample.
- **z_max**: Upper redshift bound for the population sample.

### GTI filtering
- **gti_pn_threshold**: PN count rate threshold (ct/s) above which time is flagged as a soft proton flare.
- **gti_mos_threshold**: MOS count rate threshold (ct/s) above which time is flagged as a soft proton flare.
- **gti_pn_energy_min**: Lower PI bound of the PN energy band used to build the flare-detection light curve.
- **gti_pn_energy_max**: Upper PI bound of the PN energy band used to build the flare-detection light curve.
- **gti_mos_energy_min**: Lower PI bound of the MOS energy band used to build the flare-detection light curve.
- **gti_mos_energy_max**: Upper PI bound of the MOS energy band used to build the flare-detection light curve.
- **gti_lc_binsize**: Time bin size (s) used for the flare-detection light curve.

### Extraction regions
- **src_radius**: Source extraction radius, in XMM sky pixels (1 px = 0.05").
- **bkg_radius_in**: Background annulus inner radius, in XMM sky pixels.
- **bkg_radius_out**: Background annulus outer radius, in XMM sky pixels.
- **exclude_radius**: Radius, in XMM sky pixels, excluded around any catalogued source found inside the background annulus before extraction — prevents field-specific contaminating sources from biasing the background estimate (see Known Limitations).

### Spectral & light curve processing
- **photon_index**: Assumed photon index Γ used for population-level flux/luminosity conversions.
- **channel_min_spec**: Lower limit of channel for spectroscopy.
- **channel_max_spec**: Upper limit of channel for spectroscopy.
- **bin_size_spec**: Bin size for spectroscopy.
- **energy_min_lc**: Lower energy limit for light curve.
- **energy_max_lc**: Upper energy limit for light curve.
- **bin_size_lc**: Bin size for light curve.

### Plotting
- **plot_energy_min**: Lower energy bound (keV) shown on spectral plots.
- **plot_energy_max**: Upper energy bound (keV) shown on spectral plots.
- **plot_spectra_bin_size**: Number of channels grouped together per plotted bin on spectral plots.
- **gamma_distribution_bins**: Number of bins used for the photon index distribution plot.

### Population Δαox calculation
- **xmm_band3_centre**: Assumed centre energy (keV) of 4XMM band 3, used when interpolating to the target energy.
- **xmm_band3_width**: Assumed width (keV) of 4XMM band 3, used in the same interpolation.
- **xmm_band4_centre**: Assumed centre energy (keV) of 4XMM band 4, used when interpolating to the target energy.
- **xmm_band4_width**: Assumed width (keV) of 4XMM band 4, used in the same interpolation.
- **xmm_target_energy**: Rest-frame energy (keV) Δαox is computed at (standard is 2 keV).
- **quality_sum_flag_max**: Maximum acceptable 4XMM `SC_SUM_FLAG` value; sources at or above this are excluded.
- **quality_psfflux_min**: Minimum SDSS i-band PSFFLUX (nanomaggies) required to keep a source.
- **quality_colour_cut**: Maximum allowed deviation of (g−i) colour from the sample median, used to exclude reddened/contaminated sources.
- **radio_flux_threshold**: FIRST radio flux (mJy) above which a matched source is excluded as radio-loud.
- **civ_blueshift_min**: Lower physical bound (km/s) accepted for CIV blueshift.
- **civ_blueshift_max**: Upper physical bound (km/s) accepted for CIV blueshift.
- **civ_bins**: Number of bins used for the binned CIV correlation plot.
- **del_aox_distribution_binwidth**: Fixed bin width used for the population Δαox distribution plot (kept as a fixed width rather than a fixed count so it's comparable across runs and to literature bin widths chosen to match measurement precision).
- **xray_weak_threshold**: Δαox value below which a source is classified as X-ray weak.
- **xray_strong_threshold**: Δαox value above which a source is classified as X-ray strong.

### Hardness ratio & BEHR
- **hr_soft_min**: Lower PI bound of the soft band used for hardness ratio.
- **hr_soft_max**: Upper PI bound of the soft band used for hardness ratio.
- **hr_hard_min**: Lower PI bound of the hard band used for hardness ratio.
- **hr_hard_max**: Upper PI bound of the hard band used for hardness ratio.
- **behr_nsim**: Number of Gibbs sampler simulations BEHR runs.
- **behr_nburnin**: Number of burn-in iterations discarded before BEHR starts recording samples.

### Fractional variability
- **fvar_min_bins**: Minimum number of valid light curve bins required to attempt an Fvar calculation.

### Spectral fitting
- **specfit_energy_min**: Lower energy bound (keV) used in the XSPEC fit.
- **specfit_energy_max**: Upper energy bound (keV) used in the XSPEC fit.
- **specfit_mincounts**: Minimum counts per bin used when grouping spectra for fitting.
- **specfit_gamma_guess**: Initial photon index guess given to the fit.
- **specfit_gamma_min**: Lower bound allowed on the fitted photon index (hard limit; the fit cannot converge below this).
- **specfit_gamma_max**: Upper bound allowed on the fitted photon index (hard limit; the fit cannot converge above this).
- **best_fit_min_dof**: Minimum degrees of freedom a candidate fit must have to be preferred when selecting the best-fitting epoch for a multi-ObsID target (selection itself is by minimum reduced-cstat). If no candidate clears this floor, the best available fit is used instead and logged as a fallback.

### run_steps (set each to `true` or `false`)
- **download_XMM**: Downloads datasets based on target selection rules from XMM observations.
- **download_catalogue**: Downloads the 4XMM_slim_DR14 and SDSS_DR16 catalogues.
- **pre_processing**: Runs Instrume processing, SAS and GTI filtering.
- **spectral_processing**: Runs spectral processing.
- **light_curve_processing**: Runs light curve processing.
- **population_analysis**: Runs crossmatch and calculates $\Delta\alpha_{\text{ox}}$.
- **LX_LUV_plot**: Plots UV against X-ray.
- **civ_and_hubble_plots**: Plots CIV and Hubble diagram plots.
- **spec_lc_plotting**: Plots spectral and light curve plots.
- **post_analysis**: Calculates hardness ratio, $F_{\mathrm{var}}$ and spectral fitting.
- **xray_weakness_comparison**: Runs the x-ray weakness comparison.
- **post_analysis_plots**: Generates all photon-index, hardness-ratio, and counts-based diagnostic plots from the x-ray weakness comparison output.

---

## Target Catalogue Data
Note: Catalogue currently contains minimal entries, therefore the user may need to input their own data into this for use.
- **z**: Redshift of target.
- **NH_Gal**: Milky Way hydrogen column density along the line of sight.
- **reference**: Reference(s) for the values described above.
- **notes**: Section for any general notes on target.

## Field Catalogue Data
Note: Catalogue currently contains minimal entries, therefore the user may need to input their own data into this for use.
- **H0**: Hubble constant.
- **omega_m**: Matter density parameter of the universe.
- **omega_lambda**: Dark energy (cosmological constant) density parameter of the universe.
- **reference**: Reference(s) for the values described above.

---

## Known Limitations
A full, detailed log of methodological approximations, tested-and-ruled-out hypotheses, and open issues is kept in `known_limitations.md` and is intended to be read alongside any results from this pipeline. Headline points:

- **Population Δαox** uses a fixed assumed photon index (`photon_index`) and an i-band UV proxy that drifts slightly across the sample's redshift range; several standard literature quality cuts (individual Γ filtering, UV reddening slopes, Eddington bias correction) are not implemented, giving a larger Δαox scatter (~0.33 dex) than fully-cut literature samples (~0.15 dex).
- **Background region contamination**: the fixed source-circle/background-annulus geometry can pick up field-specific contaminating soft X-ray sources that a fixed-radius annulus doesn't avoid on its own. Confirmed directly via a 4XMM crossmatch of the background region for one target. Fixed by excluding a small region around any catalogued source found inside the annulus before extraction (`_background_region`, controlled by `exclude_radius`).
- **Spectral fitting (Γ) had several real bugs, now fixed**: background was not being subtracted in the XSPEC fit at all (missing `s.background`); the reported Γ was sometimes read from a stale model state left over after a failed `Fit.error()` search rather than the true best fit; multi-instrument (combined) fits didn't propagate NH, energy band, redshift, and normalisation to the MOS model copies, leaving them at XSPEC's library defaults; spectra with negative net counts after background subtraction (an invalid detection) could still produce a boundary-pinned, physically meaningless Γ. All four are fixed; a validity gate now skips any fit with non-positive net counts rather than reporting one. Post-fix, individually-fitted Γ values now match their own confidence intervals and agree with published values (Nardini et al. 2019) for both cross-checked sources.
- **`PhoIndex` hard-limit boundary bug (found post-deployment, now fixed)**: the hard limits passed to XSPEC's `PhoIndex.values` were set half a unit beyond the intended `specfit_gamma_min`/`specfit_gamma_max` (i.e. `gamma_min − 0.5` and `gamma_max + 0.5`) rather than at those values directly. Poorly-constrained fits could converge to physically meaningless photon indices (Γ < 0, Γ > 5.5) outside the intended range, while still returning a zero-width confidence interval that superficially resembled a genuinely hard-pinned (and therefore filterable) fit. Hard limits now collapse exactly to `specfit_gamma_min`/`specfit_gamma_max`.
- **Best-fit epoch selection**: for targets with multiple archival ObsIDs, the "best" fit per instrument (PN/combined/MOS) was originally chosen by maximum degrees of freedom — i.e. the largest/most-exposed spectrum — which in practice sometimes selected the statistically worst-fitting epoch over a smaller but genuinely better-constrained one. Selection is now based on minimum reduced-cstat (`cstat/dof`), preferring candidates above `best_fit_min_dof` first and falling back to the best available fit (logged when triggered) if none clear that floor.
- **Hardness ratio `ZeroDivisionError`**: `_hardness_ratio` guarded against a zero background `BACKSCAL` but not a zero source `BACKSCAL`; either could occur for a given ObsID and crash the run. Both are now guarded, with the affected ObsID skipped and logged.
- **Photon index (Γ) can still be poorly constrained for very low-count sources** even with the above fixed — this reflects genuine Poisson-limited data, not a code issue, and is consistent with count levels the literature itself treats as marginal (Shlentsova et al. 2026 report similarly low EPIC counts, ~30–70, for some individually-studied sources). A pre-fit `_count_check` step now flags/optionally halts processing for targets below a configurable count threshold (`min_counts`) before committing to a full reduction.
- **MOS-specific Γ–HR decorrelation**: photon index correlates cleanly with hardness ratio for PN and combined fits (R²≈0.67–0.70) but not for MOS (R²≈0.003), despite MOS formal error bars being comparable to or tighter than PN's at matched counts — ruling out simple low-count noise as the explanation. Consistent with the optical-loading effect reported for MOS spectra of blue, X-ray-faint sources in Leighly et al. (2007). Not correctable within the pipeline; MOS-only Γ should be treated as less reliable than PN or combined for any HR-based conclusion.
- **Fractional variability (Fvar) is unreliable below roughly dof≈200–300** in the combined spectral fit (used here only as a rough proxy for light-curve count level): below this, Fvar values are widely scattered or recorded as upper limits, converging to a tight, low range above it. This follows directly from the excess-variance estimator's definition — at low counts, the per-bin measurement-error term subtracted from the raw variance is itself poorly determined, so the result is noise-dominated rather than reflecting true source variability.
- **Unresolved secondary excess in the population Δαox distribution** at Δαox ≈ −0.6 to −0.8, visible once the distribution plot's bin width was set to match the underlying measurement precision. Redshift, (g−i) colour, and CIV absorption-index detection rate were checked and ruled out as drivers — each is consistent with the full sample's values. Cause not yet identified.
- **Per-target result staleness**: the pipeline processes one target per run (`target_name`), and `xray_weakness_comparison` aggregates results by reading whichever `spectral_fits_<target>.csv`/`hardness_ratios_<target>.csv` files already exist on disk for every catalogued target — it does not re-run or verify them. A target can therefore silently carry results computed under an older version of the fitting code if it hasn't been individually reprocessed since a fix landed, with no automated check to catch this. Encountered directly during development. Any results table drawn from this pipeline should confirm per-target result files postdate the code version being reported.
- **PN light curves** occasionally fail during `epiclccorr` (a SAS-side segmentation fault linked to non-standard PI ranges); affected instruments are skipped automatically and `Fvar` is left uncomputed for that case.
- **Multi-target parallel runs**: `_Reprocessing()` deletes the shared `Catalogues/*.fits` catalogues unconditionally on construction when `reprocess` is truthy and `population_analysis` is in `run_steps`; the deletion isn't scoped to the triggering target, so running multiple targets in parallel with those settings can let one terminal wipe the catalogue while another is reading it. Not currently a live risk since targets are run serially — would need the catalogue reset pulled into its own explicit step to be safe under parallel runs.
- Full details, diagnostic numbers, and reasoning for all of the above are in `known_limitations.md`.

---
This program was developed by Caden Phillips, Third Year MPhys Astrophysics Student, University of Liverpool.
