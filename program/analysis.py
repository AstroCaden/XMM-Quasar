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


def _hardness_ratio(self):
    
    for obsid in self.ObsID_current:
        fits_dir = self._ObsID_Paths(obsid)["fits_dir"]

        if not (fits_dir / "pn_src_spec.fits").exists():
            if self.auto_resolve:
                self.logger.info(f"{obsid}: PN spectra not found, running spectral_processing")
                spectral_processing(self, self.channel_min_spec, self.channel_max_spec, self.bin_size_spec, self.CCF_path)

            else:
                self.logger.info(f"{obsid}: PN spectra not found, exiting")
                sys.exit(0)

    HR_results = []

    for obsid in self.ObsID_current:
        fits_dir = self._ObsID_Paths(obsid)["fits_dir"]
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


        area_ratio = src_backscal / bkg_backscal

        soft = (channel >= 500) & (channel < 2000)
        hard = (channel >= 2000) & (channel < 8000)

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

        out_path = self.star_folder / "hardness_ratios.csv"

        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=HR_results[0].keys())
            writer.writeheader()
            writer.writerows(HR_results)

        self.logger.info(f"Hardness ratios saved to {out_path}")
                




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
            if mask.sum() < 3:
                self.logger.info(f"{obsid} {name}: too few bins for Fvar")
                continue


            rate = rate[mask]
            error = error[mask]

            mu = rate.mean()
            S2 = rate.var(ddof = 1)
            sigma2_err = (error**2).mean()

            fvar_squared = (S2 - sigma2_err)/(mu**2)


            
            if fvar_squared <= 0:
                self.logger.info(f"{obsid} {name}: Fvar² negative, source not significantly variable")
                continue

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
            out_path = self.star_folder / "fvar.csv"
            with open(out_path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=results[0].keys())
                writer.writeheader()
                writer.writerows(results)
            self.logger.info(f"Fvar saved to {out_path}")
        else:
            self.logger.info("No valid Fvar measurements")


def _spectral_fitting(self):

    NH_xspec = self.NH_Gal / 1e22
    cosmo = FlatLambdaCDM(H0=self.H0, Om0=self.omega_m)
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
        xspec.AllData.clear()
        xspec.AllModels.clear()
        xspec.Xset.chatter = 0
        xspec.Xset.logChatter = 0
        xspec.Fit.statMethod = "cstat"
        xspec.Fit.query = "yes"

        
        n = 0
        spec_list = []
        for src, rmf, arf, bkg in [
            (src_pn,   fits_dir / "pn.rmf",   fits_dir / "pn.arf",   fits_dir / "pn_bkg_spec.fits"),
            (src_mos1, fits_dir / "mos1.rmf", fits_dir / "mos1.arf", fits_dir / "mos1_bkg_spec.fits"),
            (src_mos2, fits_dir / "mos2.rmf", fits_dir / "mos2.arf", fits_dir / "mos2_bkg_spec.fits"),
        ]:
            if not src.exists():
                continue
            with fits.open(src) as hdul:
                exposure = hdul["SPECTRUM"].header.get("EXPOSURE", 0)
            if exposure <= 0:
                self.logger.info(f"{obsid}: {src.name} has zero exposure, skipping")
                continue
            spec_list.append((src, rmf, arf, bkg))
            n += 1

        if n == 0:
            self.logger.info(f"{obsid}: no valid spectra with exposure, skipping fit")
            continue

        data_str = " ".join([f"{i+1}:{i+1} {src}" for i, (src, rmf, arf, bkg) in enumerate(spec_list)])
        xspec.AllData(data_str)

        for i, (src, rmf, arf, bkg) in enumerate(spec_list):
            s = xspec.AllData(i + 1)
            s.response = str(rmf)
            s.response.arf = str(arf)
            if bkg.exists():
                s.background = str(bkg)
            s.group("data mincounts 1")
            s.ignore("**-0.5 8.0-**")


        m = xspec.Model("constant*phabs*zpowerlw")
        m.phabs.nH.values = (NH_xspec, -1)
        m.zpowerlw.Redshift.values = (self.z_target, -1)
        m.zpowerlw.PhoIndex.values = (1.8, 0.1, -1.0, -1.0, 5.0, 5.0)
        m.constant.factor.values = (1, -1)
        for i in range(2, n + 1):
            xspec.AllModels(i).constant.factor.values = (1, 0.01, 0.1, 0.1, 10.0, 10.0)
            xspec.AllModels(i).phabs.nH.link = "p1"
            xspec.AllModels(i).zpowerlw.Redshift.link = "p3"
            xspec.AllModels(i).zpowerlw.PhoIndex.link = "p4"
        xspec.Fit.nIterations = 1000
        xspec.Fit.perform()
        xspec.Fit.error("1.0 4")
        gamma    = m.zpowerlw.PhoIndex.values[0]
        gamma_lo = m.zpowerlw.PhoIndex.error[0]
        gamma_hi = m.zpowerlw.PhoIndex.error[1]
        norm     = m.zpowerlw.norm.values[0]
        cstat    = xspec.Fit.statistic
        dof      = xspec.Fit.dof
        d_L    = cosmo.luminosity_distance(self.z_target).to(u.cm).value
        f_2kev = norm * (2.0) ** (1.0 - gamma) * 1e-14 / 2.418e17
        L_2kev = 4 * np.pi * d_L**2 * f_2kev * (1 + self.z_target) ** (gamma - 2)
        xspec.AllModels.calcFlux("1.0 10.0")
        flux_1_10 = xspec.AllData(1).flux[0]
        self.logger.info(f"{obsid}: Γ = {gamma:.2f} ({gamma_lo:.2f}, {gamma_hi:.2f}) cstat/dof = {cstat:.1f}/{dof}")
        results.append({
            "obsid": obsid,
            "gamma": gamma,
            "gamma_lo": gamma_lo,
            "gamma_hi": gamma_hi,
            "flux_1_10_keV": flux_1_10,
            "log_L_2kev": np.log10(L_2kev),
            "cstat": cstat,
            "dof": dof,
            "n_spectra": n,
        })
    xspec.AllData.clear()
    xspec.AllModels.clear()
    if results:
        out_path = self.star_folder / "spectral_fits.csv"
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
