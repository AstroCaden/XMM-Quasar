import numpy as np
from astropy.io import fits
import subprocess
from pathlib import Path
import sys
from program.processing import spectral_processing
import csv
from program.processing import light_curve_processing
import xspec
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import json
from astropy.coordinates import SkyCoord
from astropy.table import Table
from program.population import Population_Analysis

def _hardness_ratio(self):

    HR_results = []

    for obsid in self.ObsID_current:
        fits_dir = self._ObsID_Paths(obsid)["fits_dir"]

        if not (fits_dir / "pn_src_spec.fits").exists() or not (fits_dir / "pn_bkg_spec.fits").exists():
            if self.auto_resolve:
                self.logger.info(f"{obsid}: PN spectra not found, running spectral_processing")
                spectral_processing(self, self.channel_min_spec, self.channel_max_spec, self.bin_size_spec, self.CCF_path)
                if not (fits_dir / "pn_src_spec.fits").exists() or not (fits_dir / "pn_bkg_spec.fits").exists():
                    self.logger.info(f"{obsid}: PN spectra still missing after spectral_processing, skipping")
                    continue
            else:
                self.logger.info(f"{obsid}: PN spectra not found, exiting")
                sys.exit(0)

        src_spec = fits_dir / "pn_src_spec.fits"
        bkg_spec = fits_dir / "pn_bkg_spec.fits"
        behr_bin = self.base_dir / "BEHR"
        output_path = fits_dir / f"behr_{obsid}"

        with fits.open(src_spec) as hdul:
            data = hdul["SPECTRUM"].data
            channel = data["CHANNEL"].astype(int)
            src_counts = data["COUNTS"].astype(float)
            src_backscal = hdul["SPECTRUM"].header["BACKSCAL"]

        with fits.open(bkg_spec) as hdul:
            data = hdul["SPECTRUM"].data
            bkg_counts = data["COUNTS"].astype(float)
            bkg_backscal = hdul["SPECTRUM"].header["BACKSCAL"]

        if bkg_backscal == 0:
            self.logger.info(f"{obsid}: background BACKSCAL is zero, skipping BEHR")
            continue

        if src_backscal == 0:
            self.logger.info(f"{obsid}: source BACKSCAL is zero, skipping BEHR")
            continue

        area_ratio = bkg_backscal / src_backscal

        channel = data["CHANNEL"].astype(int)
        energy_ev = self.channel_min_spec + channel * self.bin_size_spec

        soft = (energy_ev >= self.hr_soft_min) & (energy_ev < self.hr_soft_max)
        hard = (energy_ev >= self.hr_hard_min) & (energy_ev < self.hr_hard_max)

        src_soft = int(src_counts[soft].sum())
        src_hard = int(src_counts[hard].sum())
        bkg_soft = int(bkg_counts[soft].sum())
        bkg_hard = int(bkg_counts[hard].sum())

        if src_soft + src_hard == 0:
            self.logger.info(f"{obsid}: no PN source counts, skipping BEHR")
            continue

        subprocess.run([
            str(behr_bin),
            f"softsrc={src_soft}",
            f"hardsrc={src_hard}",
            f"softbkg={bkg_soft}",
            f"hardbkg={bkg_hard}",
            f"softarea={area_ratio}",
            f"hardarea={area_ratio}",
            "algo=gibbs",
            f"nsim={self.behr_nsim}",
            f"nburnin={self.behr_nburnin}",
            "outputHR=true",
            f"output={output_path}",
        ], cwd=fits_dir, check=True)

        HR_file = Path(f"{output_path}_HR.txt")

        with open(HR_file) as f:
            lines = [l for l in f.readlines() if not l.startswith("#") and l.strip()]

        vals = lines[0].split()

        HR_results.append({
            "obsid": obsid,
            "HR_mode": float(vals[0]),
            "HR_mean": float(vals[1]),
            "HR_median": float(vals[2]),
            "HR_lo": float(vals[3]),
            "HR_hi": float(vals[4]),
            "src_soft": src_soft,
            "src_hard": src_hard,
            "bkg_soft": bkg_soft,
            "bkg_hard": bkg_hard,
        })

        self.logger.info(
            f"{obsid}: HR = {float(vals[2]):.3f} "
            f"({float(vals[3]):.3f}, {float(vals[4]):.3f})"
        )

    if HR_results:
        out_path = self.star_folder / f"hardness_ratios_{self.star_folder.name}.csv"
        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=HR_results[0].keys())
            writer.writeheader()
            writer.writerows(HR_results)
        self.logger.info(f"Hardness ratios saved to {out_path}")
    else:
        self.logger.info("No hardness ratios produced")
                




def _Fvar(self):

    results = []

    for obsid in self.ObsID_current:
        fits_dir = self._ObsID_Paths(obsid)["fits_dir"]

        instruments = ["pn", "mos1", "mos2"]

        for name in instruments:
            lc_file = fits_dir / f"{name}_corrected_lc.ds"
            if not lc_file.exists():
                if self.auto_resolve:
                    self.logger.info(f"{obsid} {name}: corrected light curve not found, running light_curve_processing")
                    light_curve_processing(self, self.energy_min_lc, self.energy_max_lc, self.bin_size_lc, self.CCF_path)
                    if not lc_file.exists():
                        self.logger.info(f"{obsid} {name}: corrected light curve still missing, skipping")
                        continue
                else:
                    self.logger.info(f"{obsid} {name}: corrected light curve not found, terminating")
                    sys.exit(0)

                


            with fits.open(lc_file) as hdul:
                data = hdul["RATE"].data
                rate = data["RATE"]
                error = data["ERROR"]

            mask = np.isfinite(rate) & np.isfinite(error) & (rate > 0)
            if mask.sum() < self.fvar_min_bins:
                self.logger.info(f"{obsid} {name}: too few bins for Fvar")
                continue


            rate = rate[mask]
            error = error[mask]

            mu = rate.mean()
            S2 = rate.var(ddof = 1)
            sigma2_err = (error**2).mean()

            fvar_squared = (S2 - sigma2_err)/(mu**2)


            
            if fvar_squared <= 0:
                self.logger.info(f"{obsid} {name}: Fvar² negative, recording as upper limit")
                fvar_val = 0.0
                fvar_err = np.sqrt(2 / mask.sum()) * sigma2_err / mu**2
            else:
                fvar_val = np.sqrt(fvar_squared)
                fvar_err = np.sqrt((np.sqrt(2 / mask.sum()) * sigma2_err / mu**2)**2 +
                                   (np.sqrt(sigma2_err / mask.sum()) * 2 * fvar_val / mu)**2)

            self.logger.info(f"{obsid} {name}: Fvar = {fvar_val:.4f} ± {fvar_err:.4f}")
            results.append({
                "obsid": obsid,
                "instrument": name,
                "Fvar": fvar_val,
                "Fvar_err": fvar_err,
            })
        if results:
            out_path = self.star_folder / f"fvar_{self.star_folder.name}.csv"
            with open(out_path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=results[0].keys())
                writer.writeheader()
                writer.writerows(results)
            self.logger.info(f"Fvar saved to {out_path}")
        else:
            self.logger.info("No valid Fvar measurements")


def _spectral_fitting(self):

    results = []

    for obsid in self.ObsID_current:
        fits_dir = self._ObsID_Paths(obsid)["fits_dir"]

        src_pn   = fits_dir / "pn_src_spec.fits"
        src_mos1 = fits_dir / "mos1_src_spec.fits"
        src_mos2 = fits_dir / "mos2_src_spec.fits"

        if not src_pn.exists() and not src_mos1.exists() and not src_mos2.exists():
            if self.auto_resolve:
                self.logger.info(f"{obsid}: no spectra found, running spectral_processing")
                spectral_processing(self, self.channel_min_spec, self.channel_max_spec, self.bin_size_spec, self.CCF_path)
                if not src_pn.exists() and not src_mos1.exists() and not src_mos2.exists():
                    self.logger.info(f"{obsid}: spectra still missing, skipping")
                    continue
            else:
                self.logger.info(f"{obsid}: no spectra found, terminating")
                sys.exit(0)

        pn_tuple   = (src_pn,   fits_dir / "pn.rmf",   fits_dir / "pn.arf",   fits_dir / "pn_bkg_spec.fits")
        mos1_tuple = (src_mos1, fits_dir / "mos1.rmf", fits_dir / "mos1.arf", fits_dir / "mos1_bkg_spec.fits")
        mos2_tuple = (src_mos2, fits_dir / "mos2.rmf", fits_dir / "mos2.arf", fits_dir / "mos2_bkg_spec.fits")

        combined_list = [t for t in [pn_tuple, mos1_tuple, mos2_tuple] if t[0].exists()]
        pn_list = [pn_tuple] if src_pn.exists() else []
        mos_list = [t for t in [mos1_tuple, mos2_tuple] if t[0].exists()]

        for fit_type, spec_list in [("combined", combined_list), ("pn", pn_list), ("mos", mos_list)]:
            if not spec_list:
                continue

            xspec.AllData.clear()
            xspec.AllModels.clear()
            xspec.Xset.chatter = 0
            xspec.Xset.logChatter = 0
            xspec.Fit.statMethod = "cstat"
            xspec.Fit.query = "yes"

            valid = []
            for src, rmf, arf, bkg in spec_list:
                if not rmf.exists() or not arf.exists():
                    self.logger.info(f"{obsid}: {rmf.name} or {arf.name} missing, excluding from fit")
                    continue
                with fits.open(src) as hdul:
                    exposure = hdul["SPECTRUM"].header.get("EXPOSURE", 0)
                if exposure > 0:
                    valid.append((src, rmf, arf, bkg))
            if not valid:
                continue



            
            data_str = " ".join([f"{i+1}:{i+1} {src}" for i, (src, rmf, arf, bkg) in enumerate(valid)])
            xspec.AllData(data_str)

            net_total = 0.0
            for src, rmf, arf, bkg in valid:
                with fits.open(src) as hdul:
                    src_data = hdul["SPECTRUM"].data
                    channel = src_data["CHANNEL"].astype(float)
                    src_counts = src_data["COUNTS"].astype(float)
                    src_backscal = hdul["SPECTRUM"].header["BACKSCAL"]

                if bkg.exists():
                    with fits.open(bkg) as hdul:
                        bkg_data = hdul["SPECTRUM"].data
                        bkg_counts = bkg_data["COUNTS"].astype(float)
                        bkg_backscal = hdul["SPECTRUM"].header["BACKSCAL"]
                    if bkg_backscal > 0:
                        area_ratio = src_backscal / bkg_backscal
                    else:
                        area_ratio = 0
                else:
                    bkg_counts = np.zeros_like(src_counts)
                    area_ratio = 0


                    

                energy_kev = (self.channel_min_spec + channel * self.bin_size_spec) / 1000.0
                band_mask = (energy_kev >= self.specfit_energy_min) & (energy_kev <= self.specfit_energy_max)
                net_total += (src_counts[band_mask] - bkg_counts[band_mask] * area_ratio).sum()


            if net_total <= 0:
                self.logger.info(f"{obsid} {fit_type}: net counts <= 0 ({net_total:.1f}) in {self.specfit_energy_min}-{self.specfit_energy_max} keV, not a valid detection, skipping fit")
                continue

            for i, (src, rmf, arf, bkg) in enumerate(valid):
                if not rmf.exists() or not arf.exists():
                    self.logger.info(f"{obsid}: {rmf.name} or {arf.name} missing, skipping this spectrum")
                    continue
                
                s = xspec.AllData(i + 1)
                s.response = str(rmf)
                s.response.arf = str(arf)
                if bkg.exists():
                    with fits.open(bkg) as hdul:
                        bkg_exposure = hdul["SPECTRUM"].header.get("EXPOSURE", 0)
                    if bkg_exposure > 0:
                        s.background = str(bkg)
                        
                        
                s.group(f"data mincounts {self.specfit_mincounts}")
                s.ignore(f"**-{self.specfit_energy_min} {self.specfit_energy_max}-**")
    
                
            NH_xspec = self.NH_Gal / 1e22
            m = xspec.Model("phabs*cflux*zpowerlw")
            m.phabs.nH.values = (NH_xspec, -1)
            m.cflux.Emin.values = (self.specfit_energy_min, -1)
            m.cflux.Emax.values = (self.specfit_energy_max, -1)
            m.zpowerlw.Redshift.values = (self.z_target, -1)
            m.zpowerlw.norm.values = (1, -1)
            
            m.zpowerlw.PhoIndex.values = (
                self.specfit_gamma_guess, 0.1,
                self.specfit_gamma_min, self.specfit_gamma_min,
                self.specfit_gamma_max, self.specfit_gamma_max,)
            
            gamma_idx = m.zpowerlw.PhoIndex.index
            
            
            flux_idx = m.cflux.lg10Flux.index
            for i in range(2, len(valid) + 1):
                mi = xspec.AllModels(i)
                mi.phabs.nH.values = (NH_xspec, -1)
                mi.cflux.Emin.values = (self.specfit_energy_min, -1)
                mi.cflux.Emax.values = (self.specfit_energy_max, -1)
                mi.zpowerlw.Redshift.values = (self.z_target, -1)
                mi.zpowerlw.norm.values = (1, -1)
                mi.zpowerlw.PhoIndex.link = f"p{gamma_idx}"

            xspec.Fit.nIterations = 1000
            try:
                xspec.Fit.perform()

                xspec.Fit.error(f"1.0 {gamma_idx}")
                gamma_lo = m.zpowerlw.PhoIndex.error[0]
                gamma_hi = m.zpowerlw.PhoIndex.error[1]

                xspec.Fit.error(f"1.0 {flux_idx}")
                log_flux_lo = m.cflux.lg10Flux.error[0]
                log_flux_hi = m.cflux.lg10Flux.error[1]

                xspec.Fit.perform()
            except Exception as e:
                self.logger.info(f"{obsid} {fit_type}: fit failed, skipping ({e})")
                continue

            gamma = m.zpowerlw.PhoIndex.values[0]
            log_flux = m.cflux.lg10Flux.values[0]
            cstat = xspec.Fit.statistic
            dof = xspec.Fit.dof

            cosmo = FlatLambdaCDM(H0=self.H0, Om0=self.omega_m)
            F = 10 ** log_flux
            emin = self.specfit_energy_min
            emax = self.specfit_energy_max
            if abs(gamma - 2.0) < 1e-6:
                K = F / (1.602176634e-9 * np.log(emax / emin))
            else:
                K = F / (1.602176634e-9 / (2 - gamma) * (emax ** (2 - gamma) - emin ** (2 - gamma)))
            f_2kev = K * (2.0) ** (1.0 - gamma) * 1.602176634e-9 / 2.418e17
            d_L = cosmo.luminosity_distance(self.z_target).to(u.cm).value
            L_2kev = 4 * np.pi * d_L**2 * f_2kev * (1 + self.z_target) ** (gamma - 2)

            self.logger.info(f"{obsid} {fit_type}: Γ = {gamma:.2f} ({gamma_lo:.2f}, {gamma_hi:.2f}) cstat/dof = {cstat:.1f}/{dof}")

            results.append({
                "obsid": obsid,
                "fit_type": fit_type,
                "gamma": gamma,
                "gamma_lo": gamma_lo,
                "gamma_hi": gamma_hi,
                "log_flux_1_10keV": log_flux,
                "log_flux_lo": log_flux_lo,
                "log_flux_hi": log_flux_hi,
                "log_L_2kev": np.log10(L_2kev),
                "cstat": cstat,
                "dof": dof,
                "counts": net_total,
                "n_spectra": len(valid),
            })

    xspec.AllData.clear()
    xspec.AllModels.clear()

    if results:
        out_path = self.star_folder / f"spectral_fits_{self.star_folder.name}.csv"
        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        self.logger.info(f"Spectral fits saved to {out_path}")
    else:
        self.logger.info("No spectral fits produced")
            

def post_analysis(self):
    _hardness_ratio(self)
    _Fvar(self)
    _spectral_fitting(self)







def _match_aox_row(self, target_name, ra, dec, aox_table):
    sdss_id = target_name.replace("SDSS_", "").replace("J", "", 1) if target_name.startswith("SDSS_J") else None

    if sdss_id:
        match = aox_table[aox_table["SDSS_NAME"] == sdss_id]
        if len(match) > 0:
            return match[0]

    coords = SkyCoord(ra=aox_table["RA_2"] * u.deg, dec=aox_table["DEC_2"] * u.deg)
    target_coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
    sep = target_coord.separation(coords)
    idx = sep.argmin()
    if sep[idx] < 5 * u.arcsec:
        return aox_table[idx]

    return None


def xray_weakness_comparison(self):
    aox_path = self.base_dir / "Catalogues" / f"aox_catalogue_{self.field_type}.fits"
    if not aox_path.exists():
        if self.auto_resolve:
            self.logger.info("aox_catalogue not found, running Population_Analysis")
            Population_Analysis(self)
            if not aox_path.exists():
                self.logger.info("aox_catalogue still missing after Population_Analysis, terminating")
                sys.exit(0)
        else:
            self.logger.info("aox_catalogue not found, terminating")
            sys.exit(0)
    with fits.open(aox_path) as hdul:
        aox_table = Table(hdul[1].data)
    target_file = self.base_dir / "target_catalogue.json"
    with open(target_file) as f:
        targets = json.load(f)
    results = []
    for name, data in targets.items():
        if data.get("type") != "science":
            continue
        try:
            pos = SkyCoord.from_name(name.replace("_", " "))
        except Exception as e:
            self.logger.info(f"{name}: name resolution failed, skipping ({e})")
            continue
 
 
        row = _match_aox_row(self, name, pos.ra.deg, pos.dec.deg, aox_table)
        if row is None:
            self.logger.info(f"{name}: no match in aox_catalogue, treating as test target")
            delta_aox = None
            classification = "test"
        else:
            delta_aox = row["DELTA_AOX"]
            if delta_aox < self.xray_weak_threshold:
                classification = "weak"
            elif delta_aox > self.xray_strong_threshold:
                classification = "strong"
            else:
                classification = "normal"
            
        safe_name = name.replace(" ", "_").replace("+", "_")
        star_folder = self.base_dir / "Datasets" / safe_name

        obsid_fits = {}
        spec_path = star_folder / f"spectral_fits_{safe_name}.csv"
        if not spec_path.exists():
            self.logger.info(f"{name}: spectral_fits_{safe_name}.csv not found, please run spectral fitting for this target first")
        else:
            with open(spec_path, newline="") as f:
                for r in csv.DictReader(f):
                    obsid_fits.setdefault(r["obsid"], {})[r["fit_type"]] = r

        obsid_hr = {}
        hr_path = star_folder / f"hardness_ratios_{safe_name}.csv"
        if not hr_path.exists():
            self.logger.info(f"{name}: hardness_ratios_{safe_name}.csv not found, please run hardness ratio calculation for this target first")
        else:
            with open(hr_path, newline="") as f:
                for r in csv.DictReader(f):
                    obsid_hr[r["obsid"]] = r

        obsid_fvar = {}
        fvar_path = star_folder / f"fvar_{safe_name}.csv"
        if not fvar_path.exists():
            self.logger.info(f"{name}: fvar_{safe_name}.csv not found, please run Fvar calculation for this target first")
        else:
            with open(fvar_path, newline="") as f:
                for r in csv.DictReader(f):
                    if r["instrument"] == "pn":
                        obsid_fvar[r["obsid"]] = r

        all_obsids = set(obsid_fits) | set(obsid_hr) | set(obsid_fvar)
        if not all_obsids:
        
            all_obsids = {None}

        def _num(row, key):
            return float(row[key]) if row is not None and row.get(key) not in (None, "") else None

        for obsid in sorted(all_obsids, key=lambda o: (o is None, o)):
            pn = obsid_fits.get(obsid, {}).get("pn")
            combined = obsid_fits.get(obsid, {}).get("combined")
            mos = obsid_fits.get(obsid, {}).get("mos")
            hr = obsid_hr.get(obsid)
            fvar = obsid_fvar.get(obsid)

            results.append({
                "target": name,
                "obsid": obsid,
                "delta_aox": delta_aox,
                "classification": classification,
                "gamma_pn": _num(pn, "gamma"),
                "gamma_pn_lo": _num(pn, "gamma_lo"),
                "gamma_pn_hi": _num(pn, "gamma_hi"),
                "gamma_combined": _num(combined, "gamma"),
                "gamma_combined_lo": _num(combined, "gamma_lo"),
                "gamma_combined_hi": _num(combined, "gamma_hi"),
                "dof_combined": _num(combined, "dof"),
                "gamma_mos": _num(mos, "gamma"),
                "gamma_mos_lo": _num(mos, "gamma_lo"),
                "gamma_mos_hi": _num(mos, "gamma_hi"),
                "dof_mos": _num(mos, "dof"),
                "hr_median": _num(hr, "HR_median"),
                "hr_lo": _num(hr, "HR_lo"),
                "hr_hi": _num(hr, "HR_hi"),
                "fvar": _num(fvar, "Fvar"),
                "fvar_err": _num(fvar, "Fvar_err"),
                "dof_pn": _num(pn, "dof"),
                "counts_pn": _num(pn, "counts"),
            
            })

        delta_aox_str = f"{delta_aox:.2f}" if delta_aox is not None else "N/A"
        self.logger.info(f"{name}: {classification} (Δαox={delta_aox_str}), {len(all_obsids)} obsid row(s) written")
 
 
    if results:
        out_path = self.base_dir / "Datasets" / f"xray_weakness_comparison_{self.field_type}.csv"
        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        self.logger.info(f"Xray weakness comparison saved to {out_path}")
    else:
        self.logger.info("No results for xray weakness comparison")
