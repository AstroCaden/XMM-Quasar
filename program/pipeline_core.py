import logging
import json
import os
import re
import shutil
import sys
from pathlib import Path
from astroquery.heasarc import Heasarc
from astropy.coordinates import SkyCoord

from program.download import Dataset_Download
from program.download import _Single_Downloader
from program.processing import pre_processing
from program.processing import spectral_processing
from program.processing import light_curve_processing
from program.plotting import plotting
from program.download import Download_Catalogue
from program.population import Population_Analysis
from program.population import lxluv_plot
from program.population import civ_and_hubble
from program.analysis import post_analysis
from program.analysis import xray_weakness_comparison

class Pipeline:

    def __init__(
        self,
        target_name,
        working_dir,
        number_of_datasets,
        obsid_to_exclude,
        obsid_specifically_chosen,
        concurrent_downloads,
        CCF_path,
        energy_min_lc,
        energy_max_lc,
        bin_size_lc,
        bin_size_spec,
        channel_min_spec,
        channel_max_spec,
        auto_resolve,
        H0,
        omega_m,
        omega_lambda,
        field_type,
        photon_index,
        civ_bins,
        z_target,
        NH_Gal,
        crossmatch_radius_arcsec,
        z_min,
        z_max,
        gti_pn_threshold,
        gti_mos_threshold,
        gti_pn_energy_min,
        gti_pn_energy_max,
        gti_mos_energy_min,
        gti_mos_energy_max,
        gti_lc_binsize,
        src_radius,
        bkg_radius_in,
        bkg_radius_out,
        hr_soft_min,
        hr_soft_max,
        hr_hard_min,
        hr_hard_max,
        behr_nsim,
        behr_nburnin,
        fvar_min_bins,
        specfit_energy_min,
        specfit_energy_max,
        specfit_mincounts,
        specfit_gamma_guess,
        specfit_gamma_min,
        specfit_gamma_max,
        specfit_const_min,
        specfit_const_max,
        quality_sum_flag_max,
        quality_psfflux_min,
        quality_colour_cut,
        civ_blueshift_min,
        civ_blueshift_max,
        plot_energy_min,
        plot_energy_max,
        plot_spectra_bin_size,
        radio_flux_threshold,
        xray_weak_threshold,
        xmm_band3_centre,
        xmm_band4_centre,
        xmm_target_energy,
        target_use,
        reprocess,
        run_steps,
        min_counts,
        stop_at_low_counts,
        exclude_radius,
        xray_strong_threshold,
    ):

        self.target_name = target_name
        self.number_of_datasets = number_of_datasets
        self.ObsID_excluded = obsid_to_exclude
        self.obsid_specifically_chosen = obsid_specifically_chosen
        self.concurrent_downloads = concurrent_downloads
        self.CCF_path = CCF_path
        self.energy_min_lc = energy_min_lc
        self.energy_max_lc = energy_max_lc
        self.bin_size_lc = bin_size_lc
        self.bin_size_spec = bin_size_spec
        self.channel_min_spec = channel_min_spec
        self.channel_max_spec = channel_max_spec
        self.auto_resolve = auto_resolve
        self.H0 = H0
        self.omega_m = omega_m
        self.omega_lambda = omega_lambda
        self.field_type = field_type
        self.photon_index = photon_index
        self.civ_bins = civ_bins
        self.z_target = z_target
        self.NH_Gal = NH_Gal
        self.crossmatch_radius_arcsec = crossmatch_radius_arcsec
        self.z_min = z_min
        self.z_max = z_max
        self.gti_pn_threshold = gti_pn_threshold
        self.gti_mos_threshold = gti_mos_threshold
        self.gti_pn_energy_min = gti_pn_energy_min
        self.gti_pn_energy_max = gti_pn_energy_max
        self.gti_mos_energy_min = gti_mos_energy_min
        self.gti_mos_energy_max = gti_mos_energy_max
        self.gti_lc_binsize = gti_lc_binsize
        self.src_radius = src_radius
        self.bkg_radius_in = bkg_radius_in
        self.bkg_radius_out = bkg_radius_out
        self.hr_soft_min = hr_soft_min
        self.hr_soft_max = hr_soft_max
        self.hr_hard_min = hr_hard_min
        self.hr_hard_max = hr_hard_max
        self.behr_nsim = behr_nsim
        self.behr_nburnin = behr_nburnin
        self.fvar_min_bins = fvar_min_bins
        self.specfit_energy_min = specfit_energy_min
        self.specfit_energy_max = specfit_energy_max
        self.specfit_mincounts = specfit_mincounts
        self.specfit_gamma_guess = specfit_gamma_guess
        self.specfit_gamma_min = specfit_gamma_min
        self.specfit_gamma_max = specfit_gamma_max
        self.specfit_const_min = specfit_const_min
        self.specfit_const_max = specfit_const_max
        self.quality_sum_flag_max = quality_sum_flag_max
        self.quality_psfflux_min = quality_psfflux_min
        self.quality_colour_cut = quality_colour_cut
        self.civ_blueshift_min = civ_blueshift_min
        self.civ_blueshift_max = civ_blueshift_max
        self.plot_energy_min = plot_energy_min
        self.plot_energy_max = plot_energy_max
        self.plot_spectra_bin_size = plot_spectra_bin_size
        self.radio_flux_threshold = radio_flux_threshold
        self.xray_weak_threshold = xray_weak_threshold
        self.xmm_band3_centre = xmm_band3_centre
        self.xmm_band4_centre = xmm_band4_centre
        self.xmm_target_energy = xmm_target_energy
        self.target_use = target_use
        self.reprocess = reprocess
        self.run_steps = run_steps
        self.min_counts = min_counts
        self.stop_at_low_counts = stop_at_low_counts
        self.exclude_radius=exclude_radius
        self.xray_strong_threshold=xray_strong_threshold

        if working_dir == "working":
            self.base_dir = Path(__file__).resolve().parent.parent
        else:
            path = Path(working_dir)
            if not path.exists():
                raise RuntimeError("Working directory does not exist")
            self.base_dir = path.resolve()

        safe_folder = re.sub(r"[^A-Za-z0-9_.-]+", "_", target_name).strip("_")

        datasets_dir = self.base_dir / "Datasets"
        datasets_dir.mkdir(exist_ok=True)

        self.star_folder = datasets_dir / safe_folder
        self.star_folder.mkdir(exist_ok=True)

        self.obsids_dir = self.star_folder / "ObsIDs"
        self.obsids_dir.mkdir(exist_ok=True)

        print(f"Inside {self.star_folder}")

        self.heasarc = Heasarc()

        self.pos = SkyCoord.from_name(self.target_name)
        self.RA = self.pos.ra.deg
        self.Dec = self.pos.dec.deg

        print(f"RA: {self.RA}, Dec: {self.Dec}")

        self.logger = logging.getLogger(self.target_name)
        self.logger.setLevel(logging.INFO)

        formatter = logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s"
        )

        self.logger.handlers.clear()

        stream = logging.StreamHandler()
        stream.setFormatter(formatter)
        self.logger.addHandler(stream)

        logfile = self.star_folder / f"{self.target_name}_pipeline.log"
        if logfile.exists():
            logfile.unlink()

        file_handler = logging.FileHandler(logfile)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

        self.failed = {}

        failed_file = self.star_folder / "failed_obsids.json"
        if failed_file.exists():
            try:
                self.failed = json.loads(failed_file.read_text())
            except Exception:
                self.failed = {}

        self._Refresh_ObsID()
        self._Reprocessing()

    def _Refresh_ObsID(self):

        if not self.obsids_dir.exists():
            self.ObsID_array = []
            self.ObsID_current = []
            return

        failed = set(self.failed.keys())

        self.ObsID_array = sorted(
            d.name
            for d in self.obsids_dir.iterdir()
            if d.is_dir()
            and d.name[:9].isdigit()
            and d.name not in failed
        )

        self.ObsID_current = self.ObsID_array.copy()

    def _Failed_ObsID(self, obsid, reason):

        message = f"{obsid} FAILED: {reason}"
        self.logger.warning(message)

        self.failed[str(obsid)] = message

        with open(self.star_folder / "failed_obsids.json", "w") as f:
            json.dump(self.failed, f, indent=2)

    def _ObsID_Paths(self, obsid):
        obs_dir = self.obsids_dir / obsid
        fits_dir = next(obs_dir.glob(f"*_{obsid}"))

        return {
            "obs_dir": obs_dir,
            "fits_dir": fits_dir,
        }



    def _Reprocessing(self):
        if self.reprocess:
            self.logger.info("Beginning reprocessing")
            steps = self.run_steps
            for obsid in self.ObsID_current:
                paths = self._ObsID_Paths(obsid)
                fits_dir = paths["fits_dir"]

                if steps.get("pre_processing", False):
                    for pattern in ["ccf.cif", "*SUM.SAS", "*ImagingEvts.ds", "*Badpixels.ds", "*_clean.ds", "*_gti.ds", "*_bkg_lc.ds"]:
                        for f in fits_dir.glob(pattern):
                            f.unlink()
                            self.logger.info(f"Removed {f}")

                if steps.get("spectral_processing", False):
                    for pattern in ["*_src_spec.fits", "*_bkg_spec.fits", "*.rmf", "*.arf"]:
                        for f in fits_dir.glob(pattern):
                            f.unlink()
                            self.logger.info(f"Removed {f}")

                if steps.get("light_curve_processing", False):
                    for pattern in ["*_src_lc.ds", "*_bkg_lc_sci.ds", "*_corrected_lc.ds"]:
                        for f in fits_dir.glob(pattern):
                            f.unlink()
                            self.logger.info(f"Removed {f}")

                if steps.get("post_analysis", False):
                    for f in fits_dir.glob("behr_*"):
                        f.unlink()
                        self.logger.info(f"Removed {f}")

            if steps.get("population_analysis", False):
                cat_dir = self.base_dir / "Catalogues"
                for name in ["matched_catalogue.fits", f"aox_catalogue_{self.field_type}.fits"]:
                    delete_file = cat_dir / name
                    if delete_file.exists():
                        delete_file.unlink()
                        self.logger.info(f"Removed {delete_file}")

            if steps.get("post_analysis", False):
                for name in ["hardness_ratios.csv", "fvar.csv", "spectral_fits.csv"]:
                    delete_file = self.star_folder / name
                    if delete_file.exists():
                        delete_file.unlink()
                        self.logger.info(f"Removed {delete_file}")








Pipeline.Dataset_Download = Dataset_Download
Pipeline._Single_Downloader = _Single_Downloader
Pipeline.pre_processing = pre_processing
Pipeline.spectral_processing = spectral_processing
Pipeline.light_curve_processing = light_curve_processing
Pipeline.plotting = plotting
Pipeline.Download_Catalogue = Download_Catalogue
Pipeline.Population_Analysis = Population_Analysis
Pipeline.lxluv_plot = lxluv_plot
Pipeline.civ_and_hubble = civ_and_hubble
Pipeline.post_analysis = post_analysis
Pipeline.xray_weakness_comparison = xray_weakness_comparison
