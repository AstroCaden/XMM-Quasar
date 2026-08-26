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
        NH_Gal
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
