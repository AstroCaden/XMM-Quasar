import numpy as np
import os
import subprocess
from astropy.time import Time
import sys
from concurrent.futures import ThreadPoolExecutor
from astroquery.esa.xmm_newton import XMMNewton
import shutil
from pathlib import Path
import requests



def _Single_Downloader(self, obsid):
    self.logger.info(f"Downloading {obsid}")

    obsid_dir = self.obsids_dir / obsid
    obsid_dir.mkdir(parents=True, exist_ok=True)

    download_script = (
        "from astroquery.esa.xmm_newton import XMMNewton\n"
        f"XMMNewton.download_data('{obsid}', level='ODF')\n"
    )

    subprocess.run(
        [
            sys.executable,
            "-c",
            download_script,
        ],
        cwd=obsid_dir,
        check=True,
    )

    archives = list(obsid_dir.glob("*.tar.gz")) + list(obsid_dir.glob("*.tar"))

    archive = archives[0]

    self.logger.info(f"Extracting {archive.name}")

    subprocess.run(
        [
            "tar",
            "-xf",
            archive.name,
            "-C",
            obsid_dir,
        ],
        cwd=obsid_dir,
        check=True,
    )

    archive.unlink()

    for inner_tar in obsid_dir.glob("*.TAR"):
        extract_dir = inner_tar.with_suffix("")
        extract_dir.mkdir(exist_ok=True)

        new_tar = extract_dir / inner_tar.name
        shutil.move(str(inner_tar), str(new_tar))

        self.logger.info(f"Extracting {new_tar.name}")

        subprocess.run(
            [
                "tar",
                "-xf",
                new_tar.name,
            ],
            cwd=extract_dir,
            check=True,
        )

        new_tar.unlink()

    self.logger.info(f"Finished {obsid}")





def Dataset_Download(self, ObsID_wanted, ObsID_excluded, ObsID_Chosen):

    self._Refresh_ObsID()

    heasarc = self.heasarc
    table = heasarc.query_object(self.target_name, mission="xmmmaster")


    today_mjd = Time.now().mjd
    table = table[(table["time"] > 0) & (table["public_date"] <= today_mjd)]
    table.sort("time")


    if "exposure" in table.colnames:
        exp_col = "exposure"
    else:
        self.logger.info("Cannot find exposure column")
        exp_col = None

    if exp_col is not None:
        before = len(table)
        table = table[table[exp_col] >= 0.0]
        removed = before - len(table)
        if removed > 0:
            self.logger.info(f"Excluded {removed} ObsIDs with exposure < 5 seconds.")

    obsids_time = np.array(table["obsid"]).astype(str)
    _, first_idx = np.unique(obsids_time, return_index=True)
    unique_obsids = obsids_time[np.sort(first_idx)]

    if ObsID_excluded is not None and len(ObsID_excluded) > 0:
        if ObsID_excluded == "current":
            ObsID_excluded = [d.name for d in self.obsids_dir.iterdir() if d.is_dir() and d.name.isdigit()]
        ObsID_excluded = set(map(str, ObsID_excluded))

        before = list(unique_obsids)
        unique_obsids = np.array([o for o in unique_obsids if o not in ObsID_excluded], dtype=str)
        removed = sorted(set(before) - set(unique_obsids))
        self.logger.info(f"Excluded {len(removed)} ObsIDs by request.")

    if ObsID_Chosen is not None and len(ObsID_Chosen) > 0:
        ObsID_Chosen = set(map(str, ObsID_Chosen))
        unique_obsids = np.array([o for o in unique_obsids if o in ObsID_Chosen], dtype=str)
        self.logger.info(f"Whitelisted {len(unique_obsids)} ObsIDs by request.")

    if len(unique_obsids) == 0:
        self.logger.info("No ObsIDs left after filtering; nothing to download.")
        return

    already_have = set(d.name for d in self.obsids_dir.iterdir() if d.is_dir() and d.name.isdigit())
    target_total = int(ObsID_wanted)
    current_total = len(already_have)
    need = max(0, target_total - current_total)

    if need == 0:
        self.logger.info(f"Already have {current_total} ObsIDs (target {target_total}); nothing to download.")
        return

    candidates = np.array([o for o in unique_obsids if o not in already_have], dtype=str)

    if len(candidates) == 0:
        self.logger.info("No new ObsIDs available to download.")
        return

    n = min(need, len(candidates))
    idx = np.linspace(0, len(candidates) - 1, n, dtype=int)
    selected_obsids = candidates[idx]

    self.logger.info(
        f"Have {current_total}, target {target_total} => need {need}. "
        f"Will download {len(selected_obsids)} ObsIDs: {selected_obsids}"
    )


    

    obsids = list(selected_obsids)




    if len(selected_obsids) == 0:
        self.logger.info("No ObsIDs found; nothing to download.")
        return

    max_workers = self.concurrent_downloads
    self.logger.info(f"Downloading {len(selected_obsids)} ObsIDs with {max_workers} workers")

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        pool.map(self._Single_Downloader, selected_obsids)



    self._Refresh_ObsID()

    if os.environ.get("PIPELINE_AFTER_DOWNLOAD") != "1":
        os.environ["PIPELINE_AFTER_DOWNLOAD"] = "1"
        os.execv(sys.executable, [sys.executable, "-m", "program.main"] + sys.argv[1:])



def Download_Catalogue(self):
    cat_dir = self.base_dir / "Catalogues"
    cat_dir.mkdir(parents=True, exist_ok=True)

    xmm_path = cat_dir / "4XMM_slim_DR14.fits.gz"
    sdss_path = cat_dir / "SDSS_DR16_quasars.fits"

    if xmm_path.exists():
        self.logger.info("4XMM-DR14 already downloaded, skipping")
    else:
        self.logger.info("Downloading 4XMM-DR14 slimline catalogue (93 MB)...")
        url = "http://xmmssc.irap.omp.eu/Catalogue/4XMM-DR14/4XMM_slim_DR14cat_v1.0.fits.gz"
        response = requests.get(url, stream=True)
        with open(xmm_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        self.logger.info("4XMM-DR14 downloaded")

    if sdss_path.exists():
        self.logger.info("SDSS DR16 already downloaded, skipping")
    else:
        self.logger.info("Downloading SDSS DR16 quasar catalogue...")
        url = "https://data.sdss.org/sas/dr16/eboss/qso/DR16Q/DR16Q_v4.fits"
        response = requests.get(url, stream=True)
        with open(sdss_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        self.logger.info("SDSS DR16 downloaded")

        
    return xmm_path, sdss_path
