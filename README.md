# X-ray Weak AGN Code
Modular pipeline created to analyse X-ray weak AGN from XMM-Newton observations.
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

### Population Δαox calculation
- **xmm_band3_centre**: Assumed centre energy (keV) of 4XMM band 3, used when interpolating to the target energy.
- **xmm_band4_centre**: Assumed centre energy (keV) of 4XMM band 4, used when interpolating to the target energy.
- **xmm_target_energy**: Rest-frame energy (keV) Δαox is computed at (standard is 2 keV).
- **quality_sum_flag_max**: Maximum acceptable 4XMM `SC_SUM_FLAG` value; sources at or above this are excluded.
- **quality_psfflux_min**: Minimum SDSS i-band PSFFLUX (nanomaggies) required to keep a source.
- **quality_colour_cut**: Maximum allowed deviation of (g−i) colour from the sample median, used to exclude reddened/contaminated sources.
- **radio_flux_threshold**: FIRST radio flux (mJy) above which a matched source is excluded as radio-loud.
- **civ_blueshift_min**: Lower physical bound (km/s) accepted for CIV blueshift.
- **civ_blueshift_max**: Upper physical bound (km/s) accepted for CIV blueshift.
- **civ_bins**: Number of bins used for the binned CIV correlation plot.
- **xray_weak_threshold**: Δαox value below which a source is classified as X-ray weak.

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
- **specfit_gamma_min**: Lower bound allowed on the fitted photon index.
- **specfit_gamma_max**: Upper bound allowed on the fitted photon index.
- **specfit_const_min**: Lower bound on the cross-normalisation constant between instruments.
- **specfit_const_max**: Upper bound on the cross-normalisation constant between instruments.

### General
- **auto_resolve**: `true` to backtrack and run missing prerequisites automatically.

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

---

## Target Catalogue Data
Note: Catalogue currently contains minimal entries, therefore the user may need to input their own data into this for use.
- **z**: Redshift of target.
- **NH_Gal**: Milky Way hydrogen column density along the line of sight.
- **reference**: Reference(s) for the values described above.

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
- **Background region contamination**: the fixed source-circle/background-annulus geometry can pick up field-specific contaminating soft X-ray counts that a fixed-radius annulus cannot avoid. This was found to bias hardness ratios (via BEHR) toward zero regardless of the true spectral shape, confirmed by a background-free diagnostic test. Any HR result accompanied by BEHR's "Monte Carlo simulations may be unstable" warning should not be treated as reliable.
- **Photon index (Γ) from individual spectral fitting** can disagree with literature values by up to ~1-2σ for low-count sources (order tens of net PN counts). A systematic investigation ruled out NH, model setup, cross-instrument linking, background subtraction, fit non-convergence, and source mis-centring as causes, and validated the core fitting pipeline against 3C 273. The discrepancy direction is source-dependent (one source too flat, one too steep), consistent with Poisson noise dominating the fit at low counts rather than a systematic pipeline error. Γ from sources with fewer than ~100 net PN counts should be treated as illustrative, not precise.
- **PN light curves** occasionally fail during `epiclccorr` (a SAS-side segmentation fault linked to non-standard PI ranges); affected instruments are skipped automatically and `Fvar` is left uncomputed for that case.
- Full details, diagnostic numbers, and reasoning for all of the above are in `known_limitations.md`.
---
This program was developed by Caden Phillips, Third Year MPhys Astrophysics Student, University of Liverpool.
