import json
import sys
import random

from program.pipeline_core import Pipeline


if __name__ == "__main__":

    config_file = sys.argv[1] if len(sys.argv) > 1 else "config.json"
    field_file = sys.argv[2] if len(sys.argv) > 2 else "field_catalogue.json"
    target_file = sys.argv[3] if len(sys.argv) > 3 else "target_catalogue.json"

    with open(config_file) as f:
        config = json.load(f)
    with open(field_file) as f:
        field = json.load(f)
    with open(target_file) as f:
        target = json.load(f)
        
    field_type = config["field"]
    field_data = field[field_type]
    
    target_type = config["target_name"]
    target_data = target[target_type]

    steps = config.get("run_steps", {})


    pipeline = Pipeline(
        target_name=config["target_name"],
        working_dir=config["working_dir"],
        number_of_datasets=config["number_of_datasets"],
        obsid_to_exclude=config["ObsID_to_exclude"],
        obsid_specifically_chosen=config["ObsID_specifically_chosen"],
        concurrent_downloads=config["concurrent_downloads"],
        CCF_path=config["CCF_path"],
        channel_min_spec=config["channel_min_spec"],
        channel_max_spec=config["channel_max_spec"],
        energy_min_lc=config["energy_min_lc"],
        energy_max_lc=config["energy_max_lc"],
        bin_size_spec=config["bin_size_spec"],
        bin_size_lc=config["bin_size_lc"],
        auto_resolve=config["auto_resolve"],
        field_type = config["field"],
        H0=field_data["H0"],
        omega_m=field_data["omega_m"],
        omega_lambda=field_data["omega_lambda"],
        photon_index = config["photon_index"],
        civ_bins=config["civ_bins"],
        z_target=target_data["z"],
        NH_Gal=target_data["NH_Gal"],
        crossmatch_radius_arcsec=config["crossmatch_radius_arcsec"],
        z_min=config["z_min"],
        z_max=config["z_max"],
        gti_pn_threshold=config["gti_pn_threshold"],
        gti_mos_threshold=config["gti_mos_threshold"],
        gti_pn_energy_min=config["gti_pn_energy_min"],
        gti_pn_energy_max=config["gti_pn_energy_max"],
        gti_mos_energy_min=config["gti_mos_energy_min"],
        gti_mos_energy_max=config["gti_mos_energy_max"],
        gti_lc_binsize=config["gti_lc_binsize"],
        src_radius=config["src_radius"],
        bkg_radius_in=config["bkg_radius_in"],
        bkg_radius_out=config["bkg_radius_out"],
        hr_soft_min=config["hr_soft_min"],
        hr_soft_max=config["hr_soft_max"],
        hr_hard_min=config["hr_hard_min"],
        hr_hard_max=config["hr_hard_max"],
        behr_nsim=config["behr_nsim"],
        behr_nburnin=config["behr_nburnin"],
        fvar_min_bins=config["fvar_min_bins"],
        specfit_energy_min=config["specfit_energy_min"],
        specfit_energy_max=config["specfit_energy_max"],
        specfit_mincounts=config["specfit_mincounts"],
        specfit_gamma_guess=config["specfit_gamma_guess"],
        specfit_gamma_min=config["specfit_gamma_min"],
        specfit_gamma_max=config["specfit_gamma_max"],
        quality_sum_flag_max=config["quality_sum_flag_max"],
        quality_psfflux_min=config["quality_psfflux_min"],
        quality_colour_cut=config["quality_colour_cut"],
        civ_blueshift_min=config["civ_blueshift_min"],
        civ_blueshift_max=config["civ_blueshift_max"],
        plot_energy_min=config["plot_energy_min"],
        plot_energy_max=config["plot_energy_max"],
        plot_spectra_bin_size=config["plot_spectra_bin_size"],
        radio_flux_threshold=config["radio_flux_threshold"],
        xray_weak_threshold=config["xray_weak_threshold"],
        xmm_band3_centre=config["xmm_band3_centre"],
        xmm_band4_centre=config["xmm_band4_centre"],
        xmm_target_energy=config["xmm_target_energy"],
        target_use=target_data["type"],
        reprocess=config["reprocess"],
        run_steps=config["run_steps"],
        min_counts=config["min_counts"],
        stop_at_low_counts=config["stop_at_low_counts"],
        exclude_radius=config["exclude_radius"],
        xray_strong_threshold=config["xray_strong_threshold"],
        xmm_band3_width=config["xmm_band3_width"],
        xmm_band4_width=config["xmm_band4_width"],
    )

    
    if steps.get("download_XMM", False):
        pipeline.Dataset_Download(
            ObsID_wanted=config["number_of_datasets"],
            ObsID_excluded=config["ObsID_to_exclude"],
            ObsID_Chosen=config["ObsID_specifically_chosen"]
            
        )
    if steps.get("download_catalogue", False):
        pipeline.Download_Catalogue()

    if steps.get("pre_processing", False):
        pipeline.pre_processing(CCF_path=config["CCF_path"],)

    if steps.get("spectral_processing", False):
        pipeline.spectral_processing(channel_min_spec=config["channel_min_spec"],
                                 channel_max_spec=config["channel_max_spec"],
                                 bin_size_spec=config["bin_size_spec"],
                                 CCF_path=config["CCF_path"]
                               
                                 
        )

    if steps.get("light_curve_processing", False):
        pipeline.light_curve_processing(energy_min_lc=config["energy_min_lc"],
                                 energy_max_lc=config["energy_max_lc"],
                                 bin_size_lc=config["bin_size_lc"],
                                 CCF_path=config["CCF_path"]

        )


    if steps.get("population_analysis", False):
        pipeline.Population_Analysis()

    if steps.get("LX_LUV_plot", False):
        pipeline.lxluv_plot(photon_index=config["photon_index"])

    
    if steps.get("civ_and_hubble_plots", False):
        pipeline.civ_and_hubble(civ_bins=config["civ_bins"], photon_index=config["photon_index"])
        
    if steps.get("spec_lc_plotting", False):
        pipeline.plotting()

    if steps.get("post_analysis", False):
        pipeline.post_analysis()

    if steps.get("xray_weakness_comparison", False):
        pipeline.xray_weakness_comparison()

    
    
