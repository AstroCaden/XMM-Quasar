# X-ray Weak AGN Code
Modular pipeline created to run NICER data through HEASoft.

HEASoft + SAS tools + BEHR must be installed and on your PATH for the program to run.

Initiate with `chmod +x run_pipeline.sh`.

Run program as `./run_pipeline.sh` in BASH.

---

## Config parameters

- **target_name**: Name of the target you wish to analyse.
- **working_dir**: Directory for the program to operate in. Use `"working"` for the directory containing the `.py` file, or provide an absolute path.
- **field**: Universe you wish to use in analysis. Options can be found in Field_Catalogue.json.
- **ObsID_to_exclude**: ObsIDs to skip. Example: `["0060010101", "0060010102"]`.
- **ObsID_specifically_chosen**: ObsIDs to use (downloads only these). Example: `["0060010101", "0060010102"]`.
- **number_of_datasets**: How many datasets wanted in total. If fewer exist, the program downloads what’s available.
- **concurrent_downloads**: Number of dataset to download concurrently. Recommended 5-10.
- **CCF_Path** Path to CCF folder.
- **channel_min_spec**: Lower limit of channel for spectroscopy.
- **channel_max_spec**: Upper limit of channel for spectroscopy.
- **bin_size_spec**: Bin size for spectroscopy.
- **energy_min_lc**: Lower energy limit for light curve.
- **energy_max_lc**: Upper energy limit for light curve.
- **bin_size_lc**: Bin size for light curve.
- **auto_resolve**: `true` to backtrack and run missing prerequisites automatically.
- **civ_bins**: Bin size for CIV plot.

### run_steps (set each to `true` or `false`)

- **download_XMM**: Downloads datasets based on target selection rules from XMM observations.
- **download_catalogue**: Downloads the 4XMM_slim_DR14 and SDSS_DR16 catalogues.
- **pre_processing**: Runs Instrume processing, SAS and gti filtering.
- **spectral_processing**: Runs spectral processing.
- **light_curve_processing**: Runs light curve processing.
- **population_analysis**: Runs crossmatch and calculates  \(\Delta\alpha_{\text{ox}}\).
- **LX_LUV_plot**: Plots UV against X-ray.
- **civ_and_hubble_plots**: Plots civ and hubble plots.
- **spec_lc_plotting**: Plots spectrographic plot and light curve.
- **post_analysis**: Calculates hardness ratio, \(F_{\text{var}}\) and spectral fitting.
---

## Catalogue Data
Note: Catalogue currently contains minimal enteries, therefore the user may need to input their own data into this for use. 
- **orbital_period**: Orbital period of object (in days).
- **orbital_pdot**: First time derivative of orbital period (in s/s). Use `0.0` if not wanted.
- **reference_epoch**: Phase-zero epoch (**T0**) in MJD (e.g. Tasc / inferior conjunction from literature).
- **reference**: Reference(s) for the values described above.

---

This program was developed by Caden Phillips, Third Year MPhys Astrophysics Student, University of Liverpool.
