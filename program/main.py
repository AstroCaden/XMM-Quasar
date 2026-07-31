import json
import sys

from program.pipeline_core import Pipeline


if __name__ == "__main__":

    config_file = sys.argv[1] if len(sys.argv) > 1 else "config.json"
    field_file    = sys.argv[2] if len(sys.argv) > 2 else "field_catalogue.json"

    with open(config_file) as f:
        config = json.load(f)
    with open(field_file) as f:
        field = json.load(f)
        
    field_type = config["field"]
    field_data = field[field_type]

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
        H0 = field_data["H0"],
        omega_m = field_data["omega_m"],
        omega_lambda = field_data["omega_lambda"],
        field_type = config["field"],
        photon_index = config["photon_index"],
        civ_bins=config["civ_bins"]
        
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
