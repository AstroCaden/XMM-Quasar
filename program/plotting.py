from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from program.processing import light_curve_processing, spectral_processing
import csv
import sys
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import json
from program.analysis import _match_aox_row



def _plot_light_curves(self):
    for obsid in self.ObsID_current:
        fits_dir = self._ObsID_Paths(obsid)["fits_dir"]

        src_lc_files = list(fits_dir.glob("*_corrected_lc.ds"))
        if not src_lc_files:
            if self.auto_resolve:
                self.logger.info(f"{obsid}: source light curves not found, running light_curve_processing")
                light_curve_processing(self, self.energy_min_lc, self.energy_max_lc, self.bin_size_lc, self.CCF_path)
                src_lc_files = list(fits_dir.glob("*_corrected_lc.ds"))
                if not src_lc_files:
                    self.logger.info(f"{obsid}: source light curves still missing after light_curve_processing, skipping")
                    continue
            else:
                self.logger.info(f"{obsid}: source light curves not found, terminating")
                sys.exit(0)

        plots_dir = fits_dir / "Plots"
        plots_dir.mkdir(parents=True, exist_ok=True)

        instruments = [
            ("pn",   "EPIC-PN"),
            ("mos1", "EPIC-MOS1"),
            ("mos2", "EPIC-MOS2"),
        ]
        fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

        for ax, (name, label) in zip(axes, instruments):
            lc_file = fits_dir / f"{name}_corrected_lc.ds"
            if not lc_file.exists():
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
                ax.set_ylabel("Count rate (ct/s)", fontsize=11)
                ax.text(0.02, 0.95, label, transform=ax.transAxes, fontsize=11, va="top")
                continue

            with fits.open(lc_file) as hdul:
                data = hdul["RATE"].data
                time = data["TIME"]
                rate = data["RATE"]
                error = data["ERROR"]

            time = time - time[0]
            mask = np.isfinite(rate) & np.isfinite(error)

            ax.step(time[mask], rate[mask],
                    where="mid", linewidth=0.8, color="black")

            for x, y, e in zip(time[mask], rate[mask], error[mask]):
                ax.plot([x, x], [y - e, y + e],
                        color="black", linewidth=0.6, linestyle="--", alpha=0.5)
                ax.plot([x - 0.01*x, x + 0.01*x], [y - e, y - e],
                        color="black", linewidth=0.6, alpha=0.5)
                ax.plot([x - 0.01*x, x + 0.01*x], [y + e, y + e],
                        color="black", linewidth=0.6, alpha=0.5)
                ax.plot([x], [y], "o", markersize=3, color="red")

            ax.set_ylabel("Count rate (ct/s)", fontsize=11)
            ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
            ax.grid(True, which="both", alpha=0.2, linestyle="--")
            ax.text(0.02, 0.95, label, transform=ax.transAxes,
                    fontsize=11, va="top")

        axes[-1].set_xlabel("Time (s)", fontsize=12)
        fig.suptitle(f"{self.target_name}  |  ObsID {obsid}", fontsize=12)
        plt.tight_layout()
        outpath = plots_dir / f"light_curves_{obsid}.png"
        plt.savefig(outpath, dpi=150, bbox_inches="tight")
        plt.close()
        self.logger.info("Light curve plotted.")

        rows = []
        for name, label in instruments:
            lc_file = fits_dir / f"{name}_corrected_lc.ds"
            if not lc_file.exists():
                continue
            with fits.open(lc_file) as hdul:
                data = hdul["RATE"].data
                time = data["TIME"]
                rate = data["RATE"]
                error = data["ERROR"]
            time = time - time[0]
            mask = np.isfinite(rate) & np.isfinite(error)
            for t, r, err in zip(time[mask], rate[mask], error[mask]):
                rows.append([name, t, r, err])

        with open(plots_dir / f"light_curves_{obsid}.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["instrument", "time_s", "count_rate", "error"])
            writer.writerows(rows)


def _plot_spectra(self):
    for obsid in self.ObsID_current:
        fits_dir = self._ObsID_Paths(obsid)["fits_dir"]

        src_spec_files = list(fits_dir.glob("*_src_spec.fits"))
        if not src_spec_files:
            if self.auto_resolve:
                self.logger.info(f"{obsid}: source spectra not found, running spectral_processing")
                spectral_processing(self, self.channel_min_spec, self.channel_max_spec, self.bin_size_spec, self.CCF_path)
                src_spec_files = list(fits_dir.glob("*_src_spec.fits"))
                if not src_spec_files:
                    self.logger.info(f"{obsid}: source spectra still missing after spectral_processing, skipping")
                    continue
            else:
                self.logger.info(f"{obsid}: source spectra not found, terminating")
                sys.exit(0)

        plots_dir = fits_dir / "Plots"
        plots_dir.mkdir(parents=True, exist_ok=True)

        instruments = [
            ("pn",   "EPIC-PN"),
            ("mos1", "EPIC-MOS1"),
            ("mos2", "EPIC-MOS2"),
        ]

        fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

        for ax, (name, label) in zip(axes, instruments):
            spec_file = fits_dir / f"{name}_src_spec.fits"
            if not spec_file.exists():
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
                ax.set_ylabel("Counts", fontsize=11)
                ax.text(0.02, 0.95, label, transform=ax.transAxes, fontsize=11, va="top")
                continue
            with fits.open(spec_file) as hdul:
                data = hdul["SPECTRUM"].data
                channel = data["CHANNEL"].astype(float)
                counts  = data["COUNTS"].astype(float)

            energy_kev = channel / 1000.0

            bin_size = self.bin_size_spec
            n_bins = len(counts) // bin_size
            energy_binned = energy_kev[:n_bins * bin_size].reshape(n_bins, bin_size).mean(axis=1)
            counts_binned = counts[:n_bins * bin_size].reshape(n_bins, bin_size).sum(axis=1)
            errors_binned = np.sqrt(np.maximum(counts_binned, 1))

            mask = (counts_binned > 0) & (energy_binned >= self.plot_energy_min) & (energy_binned <= self.plot_energy_max)

            ax.step(energy_binned[mask], counts_binned[mask],
                    where="mid", linewidth=0.8, color="black")
            for x, y, e in zip(energy_binned[mask], counts_binned[mask], errors_binned[mask]):
                ax.plot([x, x], [y - e, y + e],
                        color="black", linewidth=0.6, linestyle="--", alpha=0.5)
                ax.plot([x - 0.01*x, x + 0.01*x], [y - e, y - e],
                        color="black", linewidth=0.6, alpha=0.5)
                ax.plot([x - 0.01*x, x + 0.01*x], [y + e, y + e],
                        color="black", linewidth=0.6, alpha=0.5)
                ax.plot([x], [y], "o", markersize=3, color="red")

            ax.set_xscale("log")
            ax.set_ylabel("Counts", fontsize=11)
            ax.set_xlim(self.plot_energy_min, self.plot_energy_max)
            ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
            ax.grid(True, which="both", alpha=0.2, linestyle="--")
            ax.text(0.02, 0.95, label, transform=ax.transAxes,
                    fontsize=11, va="top")

        axes[-1].set_xlabel("Energy (keV)", fontsize=12)
        fig.suptitle(f"{self.target_name}  |  ObsID {obsid}", fontsize=12)
        plt.tight_layout()
        outpath = plots_dir / f"spectra_{obsid}.png"
        plt.savefig(outpath, dpi=150, bbox_inches="tight")
        plt.close()
        self.logger.info("Spectra plotted.")

        all_rows = []
        for name, label in instruments:
            spec_file = fits_dir / f"{name}_src_spec.fits"
            if not spec_file.exists():
                continue
            with fits.open(spec_file) as hdul:
                data = hdul["SPECTRUM"].data
                channel = data["CHANNEL"].astype(float)
                counts = data["COUNTS"].astype(float)
            energy_kev = channel / 1000.0
            errors = np.sqrt(np.maximum(counts, 1))
            for e, c, err in zip(energy_kev, counts, errors):
                all_rows.append([name, e, c, err])

        with open(plots_dir / f"spectra_{obsid}.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["instrument", "energy_kev", "counts", "error"])
            writer.writerows(all_rows)



def _plot_combined_light_curve(self):
    plots_dir = self.star_folder / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    instruments = [
        ("pn",   "EPIC-PN"),
        ("mos1", "EPIC-MOS1"),
        ("mos2", "EPIC-MOS2"),
    ]

    all_rows = []
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=False)

    for obsid in self.ObsID_current:
        fits_dir = self._ObsID_Paths(obsid)["fits_dir"]
        lc_csv = fits_dir / "Plots" / f"light_curves_{obsid}.csv"
        if not lc_csv.exists():
            continue

        with open(lc_csv, newline="") as f:
            rows = list(csv.DictReader(f))

        for row in rows:
            all_rows.append({"obsid": obsid, **row})

        for ax, (name, label) in zip(axes, instruments):
            obs_rows = [r for r in rows if r["instrument"] == name]
            if not obs_rows:
                continue

            time = np.array([float(r["time_s"]) for r in obs_rows])
            rate = np.array([float(r["count_rate"]) for r in obs_rows])
            error = np.array([float(r["error"]) for r in obs_rows])

            ax.step(time, rate, where="mid", linewidth=0.8, color="black")
            for x, y, e in zip(time, rate, error):
                ax.plot([x, x], [y - e, y + e], color="black", linewidth=0.6, linestyle="--", alpha=0.5)
                ax.plot([x - 0.01*x, x + 0.01*x], [y - e, y - e], color="black", linewidth=0.6, alpha=0.5)
                ax.plot([x - 0.01*x, x + 0.01*x], [y + e, y + e], color="black", linewidth=0.6, alpha=0.5)
                ax.plot([x], [y], "o", markersize=3, color="red")

            ax.axvline(time[-1], color="gray", linewidth=0.5, linestyle=":")
            ax.text(time[len(time)//2], ax.get_ylim()[1] * 0.9 if ax.get_ylim()[1] != 1.0 else 0.9,
                    obsid, fontsize=7, ha="center", color="gray")

    for ax, (name, label) in zip(axes, instruments):
        ax.set_ylabel("Count rate (ct/s)", fontsize=11)
        ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        ax.grid(True, which="both", alpha=0.2, linestyle="--")
        ax.text(0.02, 0.95, label, transform=ax.transAxes, fontsize=11, va="top")

    axes[-1].set_xlabel("Time (s)", fontsize=12)
    fig.suptitle(f"{self.target_name} — combined light curves", fontsize=12)
    plt.tight_layout()
    plt.savefig(plots_dir / "combined_light_curves.png", dpi=150, bbox_inches="tight")
    plt.close()
    self.logger.info("Combined light curve plotted.")

    with open(plots_dir / f"{self.target_name}_combined_light_curves.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["obsid", "instrument", "time_s", "count_rate", "error"])
        writer.writeheader()
        writer.writerows(all_rows)

        

def _plot_combined_spectra(self):
    plots_dir = self.star_folder / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    instruments = [
        ("pn",   "EPIC-PN"),
        ("mos1", "EPIC-MOS1"),
        ("mos2", "EPIC-MOS2"),
    ]

    all_rows = []
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    for obsid in self.ObsID_current:
        fits_dir = self._ObsID_Paths(obsid)["fits_dir"]
        spec_csv = fits_dir / "Plots" / f"spectra_{obsid}.csv"
        if not spec_csv.exists():
            continue

        with open(spec_csv, newline="") as f:
            rows = list(csv.DictReader(f))

        for row in rows:
            all_rows.append({"obsid": obsid, **row})

    for ax, (name, label) in zip(axes, instruments):
        name_rows = [r for r in all_rows if r["instrument"] == name]
        if not name_rows:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
            ax.set_ylabel("Counts", fontsize=11)
            ax.text(0.02, 0.95, label, transform=ax.transAxes, fontsize=11, va="top")
            continue

        energy = np.array([float(r["energy_kev"]) for r in name_rows])
        counts = np.array([float(r["counts"]) for r in name_rows])

        bin_size = self.plot_spectra_bin_size
        n_bins = len(counts) // bin_size
        energy_binned = energy[:n_bins * bin_size].reshape(n_bins, bin_size).mean(axis=1)
        counts_binned = counts[:n_bins * bin_size].reshape(n_bins, bin_size).sum(axis=1)
        errors_binned = np.sqrt(np.maximum(counts_binned, 1))

        mask = (counts_binned > 0) & (energy_binned >= self.plot_energy_min) & (energy_binned <= self.plot_energy_max)

        ax.step(energy_binned[mask], counts_binned[mask],
                where="mid", linewidth=0.8, color="black")
        for x, y, e in zip(energy_binned[mask], counts_binned[mask], errors_binned[mask]):
            ax.plot([x, x], [y - e, y + e],
                    color="black", linewidth=0.6, linestyle="--", alpha=0.5)
            ax.plot([x - 0.01*x, x + 0.01*x], [y - e, y - e],
                    color="black", linewidth=0.6, alpha=0.5)
            ax.plot([x - 0.01*x, x + 0.01*x], [y + e, y + e],
                    color="black", linewidth=0.6, alpha=0.5)
            ax.plot([x], [y], "o", markersize=3, color="red")

        ax.set_xscale("log")
        ax.set_ylabel("Counts", fontsize=11)
        ax.set_xlim(self.plot_energy_min, self.plot_energy_max)
        ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        ax.grid(True, which="both", alpha=0.2, linestyle="--")
        ax.text(0.02, 0.95, label, transform=ax.transAxes, fontsize=11, va="top")

    axes[-1].set_xlabel("Energy (keV)", fontsize=12)
    fig.suptitle(f"{self.target_name} — combined spectra", fontsize=12)
    plt.tight_layout()
    plt.savefig(plots_dir / "combined_spectra.png", dpi=150, bbox_inches="tight")
    plt.close()
    self.logger.info("Combined spectra plotted.")

    with open(plots_dir / f"{self.target_name}_combined_spectra.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["obsid", "instrument", "energy_kev", "counts", "error"])
        writer.writeheader()
        writer.writerows(all_rows)
                    

def plotting(self):
    _plot_light_curves(self)
    _plot_spectra(self)
    _plot_combined_light_curve(self)
    _plot_combined_spectra(self)



def _plot_gamma_daox(self):

    csv_path = self.base_dir / "Datasets" / f"xray_weakness_comparison_{self.field_type}.csv"
    if not csv_path.exists():
        self.logger.info("xray_weakness_comparison csv not found, skipping Γ vs Δαox plot")
        return

    gamma_types = ["pn", "combined", "mos"]
    markers = {"pn": "x", "combined": "s", "mos": "^"}
    colors = {"weak": "blue", "normal": "black", "strong": "red"}

    data = {(cls, t): {"daox": [], "gamma": [], "lo": [], "hi": []} for cls in colors for t in gamma_types}

    with open(csv_path, newline="") as f:
        for r in csv.DictReader(f):
            if r["classification"] not in colors:
                continue
            for t in gamma_types:
                gamma_col, gamma_lo_col, gamma_hi_col = f"gamma_{t}", f"gamma_{t}_lo", f"gamma_{t}_hi"
                if not r[gamma_col]:
                    continue
                gamma_lo_val = float(r[gamma_lo_col])
                gamma_hi_val = float(r[gamma_hi_col])
                if gamma_lo_val == 0 and gamma_hi_val == 0:
                    continue
                gamma = float(r[gamma_col])
                if gamma_lo_val > gamma or gamma_hi_val < gamma:
                    self.logger.info(f"{r['target']} {t}: inverted error bounds, skipping")
                    continue
                g = data[(r["classification"], t)]
                g["daox"].append(float(r["delta_aox"]))
                g["gamma"].append(gamma)
                g["lo"].append(gamma - gamma_lo_val)
                g["hi"].append(gamma_hi_val - gamma)

    fig, ax = plt.subplots(figsize=(8, 6))
    for cls in colors:
        for t in gamma_types:
            g = data[(cls, t)]
            if not g["gamma"]:
                continue
            ax.errorbar(g["daox"], g["gamma"], yerr=[g["lo"], g["hi"]],
                         fmt=markers[t], color=colors[cls], capsize=3, markersize=6)

    fit_x = []
    fit_y = []
    for cls in colors:
        for t in gamma_types:
            fit_x.extend(data[(cls, t)]["daox"])
            fit_y.extend(data[(cls, t)]["gamma"])

    if len(fit_x) >= 2:
        slope, intercept = np.polyfit(fit_x, fit_y, 1)
        x_fit = np.linspace(min(fit_x), max(fit_x), 200)
      
        ax.plot(x_fit, slope * x_fit + intercept, color="gray", linewidth=0.8, linestyle="--", alpha=0.5)


        r2 = np.corrcoef(fit_x, fit_y)[0, 1]**2
        ax.text(0.05, 0.05, rf"$\Gamma = {slope:.3f}\Delta\alpha_{{ox}} {intercept:+.3f}$, $R^2={r2:.3f}$",
                transform=ax.transAxes, fontsize=10)

    ax.axvline(self.xray_weak_threshold, color="blue", linewidth=0.8, linestyle="--")
    ax.axvline(self.xray_strong_threshold, color="red", linewidth=0.8, linestyle="--")
    ax.set_xlabel(r"$\Delta\alpha_{ox}$", fontsize=12)
    ax.set_ylabel(r"Photon index $\Gamma$", fontsize=12)
    ax.set_title(f"Γ vs Δαox by classification | {self.field_type}", fontsize=12)

    class_handles = [Line2D([0], [0], marker="o", color=colors[cls], linestyle="None", markersize=6, label=cls) for cls in colors]
    type_handles = [Line2D([0], [0], marker=markers[t], color="gray", linestyle="None", markersize=6, label=t) for t in gamma_types]
    class_legend = ax.legend(handles=class_handles, fontsize=10, loc="upper left")
    ax.add_artist(class_legend)
    ax.legend(handles=type_handles, fontsize=10, loc="upper right")
    ax.grid(True, alpha=0.2, linestyle="--")

    plots_dir = self.base_dir / "Catalogues" / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    outpath = plots_dir / f"gamma_vs_daox_{self.field_type}.png"
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    self.logger.info(f"Γ vs Δαox plot saved to {outpath}")


def _plot_gamma_counts(self):
    
    csv_path = self.base_dir / "Datasets" / f"xray_weakness_comparison_{self.field_type}.csv"
    if not csv_path.exists():
        self.logger.info("xray_weakness_comparison csv not found, skipping Γ vs counts plot")
        return

    gamma_types = ["pn", "combined", "mos"]
    markers = {"pn": "x", "combined": "s", "mos": "^"}
    status_colors = {"clean": "black", "broken": "red"}

    data = {(status, t): {"dof": [], "gamma": []} for status in status_colors for t in gamma_types}

    with open(csv_path, newline="") as f:
        for r in csv.DictReader(f):
            for t in gamma_types:
                gamma_col, gamma_lo_col, gamma_hi_col, dof_col = f"gamma_{t}", f"gamma_{t}_lo", f"gamma_{t}_hi", f"dof_{t}"
                if not r[gamma_col] or not r[dof_col]:
                    continue
                gamma = float(r[gamma_col])
                gamma_lo = float(r[gamma_lo_col])
                gamma_hi = float(r[gamma_hi_col])
                dof = float(r[dof_col])
                status = "clean" if not (gamma_lo == 0 and gamma_hi == 0) else "broken"
                g = data[(status, t)]
                g["dof"].append(dof)
                g["gamma"].append(gamma)

    fig, ax = plt.subplots(figsize=(8, 6))
    for status in status_colors:
        for t in gamma_types:
            g = data[(status, t)]
            if not g["gamma"]:
                continue
            ax.scatter(g["dof"], g["gamma"], color=status_colors[status], marker=markers[t], s=40)

    ax.set_xscale("log")
    ax.set_xlabel("Fit degrees of freedom (proxy for counts)", fontsize=12)
    ax.set_ylabel(r"Photon index $\Gamma$", fontsize=12)
    ax.set_title(f"Fit reliability vs counts | {self.field_type}", fontsize=12)

    status_handles = [Line2D([0], [0], marker="o", color=status_colors[s], linestyle="None", markersize=6,
                              label="Well-constrained" if s == "clean" else "Boundary-pinned / broken") for s in status_colors]
    type_handles = [Line2D([0], [0], marker=markers[t], color="gray", linestyle="None", markersize=6, label=t) for t in gamma_types]
    status_legend = ax.legend(handles=status_handles, fontsize=10, loc="upper left")
    ax.add_artist(status_legend)
    ax.legend(handles=type_handles, fontsize=10, loc="upper right")
    ax.grid(True, alpha=0.2, linestyle="--")

    plots_dir = self.base_dir / "Catalogues" / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    outpath = plots_dir / f"gamma_vs_counts_{self.field_type}.png"
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    self.logger.info(f"Γ vs counts plot saved to {outpath}")




def _plot_targets_on_lxluv(self):

    aox_path = self.base_dir / "Catalogues" / f"aox_catalogue_{self.field_type}.fits"
    if not aox_path.exists():
        self.logger.info("aox_catalogue not found, skipping overlay plot")
        return

    with fits.open(aox_path) as hdul:
        aox_table = Table(hdul[1].data)

    valid = aox_table["VALID_AOX"].data.astype(bool)
    tv = aox_table[valid]
    log_L2500 = np.log10(tv["L2500"].data)
    log_L2kev = np.log10(tv["L2KEV"].data)

    target_file = self.base_dir / "target_catalogue.json"
    with open(target_file) as f:
        targets = json.load(f)

    overlay_x, overlay_y, overlay_names = [], [], []
    for name, data in targets.items():
        if data.get("type") != "science":
            continue
        try:
            pos = SkyCoord.from_name(name.replace("_", " "))
        except Exception as e:
            self.logger.info(f"{name}: name resolution failed, skipping ({e})")
            continue
        row = _match_aox_row(self, name, pos.ra.deg, pos.dec.deg, aox_table)
        if row is None or not row["VALID_AOX"]:
            continue
        overlay_x.append(np.log10(row["L2500"]))
        overlay_y.append(np.log10(row["L2KEV"]))
        overlay_names.append(name)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(log_L2500, log_L2kev, s=2, color="lightgray", alpha=0.4, label="Population")
    ax.scatter(overlay_x, overlay_y, s=25, color="red", zorder=5, label="Individually processed targets")

    ax.set_xlabel(r"$\log(L_{2500\AA})$ [erg s$^{-1}$ Hz$^{-1}$]", fontsize=12)
    ax.set_ylabel(r"$\log(L_{2\,\mathrm{keV}})$ [erg s$^{-1}$ Hz$^{-1}$]", fontsize=12)
    ax.set_title(f"Processed targets on the population LX-LUV plane | {self.field_type}", fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.2, linestyle="--")

    plots_dir = self.base_dir / "Catalogues" / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    outpath = plots_dir / f"targets_on_lxluv_{self.field_type}.png"
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    self.logger.info(f"Targets-on-LXLUV plot saved to {outpath}")






def _plot_hr_gamma(self):
 
    csv_path = self.base_dir / "Datasets" / f"xray_weakness_comparison_{self.field_type}.csv"
    if not csv_path.exists():
        self.logger.info("xray_weakness_comparison csv not found, skipping Γ vs HR plot")
        return
 
    gamma_types = ["pn", "combined", "mos"]
    markers = {"pn": "x", "combined": "s", "mos": "^"}
    colors = {"weak": "blue", "normal": "black", "strong": "red"}
 
    data = {(cls, t): {"hr": [], "hr_lo": [], "hr_hi": [], "gamma": [], "lo": [], "hi": []} for cls in colors for t in gamma_types}
    all_hr, all_gamma = [], []
 
    with open(csv_path, newline="") as f:
        for r in csv.DictReader(f):
            if r["classification"] not in colors:
                continue
            if not r["hr_median"]:
                continue
            hr = float(r["hr_median"])
            hr_lo_val = float(r["hr_lo"]) if r["hr_lo"] else hr
            hr_hi_val = float(r["hr_hi"]) if r["hr_hi"] else hr
            for t in gamma_types:
                gamma_col, gamma_lo_col, gamma_hi_col = f"gamma_{t}", f"gamma_{t}_lo", f"gamma_{t}_hi"
                if not r[gamma_col]:
                    continue
                gamma_lo_val = float(r[gamma_lo_col])
                gamma_hi_val = float(r[gamma_hi_col])
                if gamma_lo_val == 0 and gamma_hi_val == 0:
                    continue
                gamma = float(r[gamma_col])
                if gamma_lo_val > gamma or gamma_hi_val < gamma:
                    self.logger.info(f"{r['target']} {t}: inverted error bounds, skipping")
                    continue
                g = data[(r["classification"], t)]
                g["hr"].append(hr)
                g["hr_lo"].append(hr - hr_lo_val)
                g["hr_hi"].append(hr_hi_val - hr)
                g["gamma"].append(gamma)
                g["lo"].append(gamma - gamma_lo_val)
                g["hi"].append(gamma_hi_val - gamma)
                all_hr.append(hr)
                all_gamma.append(gamma)
 
    fig, ax = plt.subplots(figsize=(8, 6))
    for cls in colors:
        for t in gamma_types:
            g = data[(cls, t)]
            if not g["gamma"]:
                continue
            ax.errorbar(g["hr"], g["gamma"], xerr=[g["hr_lo"], g["hr_hi"]], yerr=[g["lo"], g["hi"]],
                         fmt=markers[t], color=colors[cls], capsize=3, markersize=6)
 
    if len(all_hr) > 1:
        slope, intercept = np.polyfit(all_hr, all_gamma, 1)
        pred = slope * np.array(all_hr) + intercept
        ss_res = np.sum((np.array(all_gamma) - pred) ** 2)
        ss_tot = np.sum((np.array(all_gamma) - np.mean(all_gamma)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot else 0
        ax.text(0.95, 0.05, f"$R^2$ = {r2:.3f}", transform=ax.transAxes, fontsize=11, ha="right", va="bottom")
 
    ax.set_xlabel("Hardness Ratio", fontsize=12)
    ax.set_ylabel(r"Photon index $\Gamma$", fontsize=12)
    ax.set_title(f"Γ vs HR by classification | {self.field_type}", fontsize=12)
 
    class_handles = [Line2D([0], [0], marker="o", color=colors[cls], linestyle="None", markersize=6, label=cls) for cls in colors]
    type_handles = [Line2D([0], [0], marker=markers[t], color="gray", linestyle="None", markersize=6, label=t) for t in gamma_types]
    class_legend = ax.legend(handles=class_handles, fontsize=10, loc="upper left")
    ax.add_artist(class_legend)
    ax.legend(handles=type_handles, fontsize=10, loc="upper right")
    ax.grid(True, alpha=0.2, linestyle="--")
 
    plots_dir = self.base_dir / "Catalogues" / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    outpath = plots_dir / f"gamma_vs_hr_{self.field_type}.png"
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    self.logger.info(f"Γ vs HR plot saved to {outpath}")




def _plot_gamma_distribution(self):

    csv_path = self.base_dir / "Datasets" / f"xray_weakness_comparison_{self.field_type}.csv"
    if not csv_path.exists():
        self.logger.info("xray_weakness_comparison csv not found, skipping Γ distribution plot")
        return

    gamma_types = ["pn", "combined", "mos"]
    colors = {"weak": "blue", "normal": "black", "strong": "red"}

    data = {cls: [] for cls in colors}

    with open(csv_path, newline="") as f:
        for r in csv.DictReader(f):
            if r["classification"] not in colors:
                continue
            for t in gamma_types:
                gamma_col, gamma_lo_col, gamma_hi_col = f"gamma_{t}", f"gamma_{t}_lo", f"gamma_{t}_hi"
                if not r[gamma_col]:
                    continue
                gamma_lo_val = float(r[gamma_lo_col])
                gamma_hi_val = float(r[gamma_hi_col])
                if gamma_lo_val == 0 and gamma_hi_val == 0:
                    continue
                gamma = float(r[gamma_col])
                if gamma_lo_val > gamma or gamma_hi_val < gamma:
                    continue
                data[r["classification"]].append((gamma, (gamma_hi_val - gamma_lo_val) / 2))

    fig, ax = plt.subplots(figsize=(8, 6))
    bins = np.linspace(0, 4, self.gamma_distrib_bins+1)
    grid = np.linspace(0, 4, 400)

    for cls in colors:
        vals = [g for g, s in data[cls]]
        if not vals:
            continue
        ax.hist(vals, bins=bins, density=True, histtype="stepfilled",
                color=colors[cls], alpha=0.25, label=f"{cls} (n={len(vals)})")

        pdf = np.zeros_like(grid)
        for g, s in data[cls]:
            if s <= 0:
                continue
            pdf += np.exp(-0.5 * ((grid - g) / s) ** 2) / (s * np.sqrt(2 * np.pi))
        if pdf.sum() > 0:
            pdf /= np.trapz(pdf, grid)
            ax.plot(grid, pdf, color=colors[cls], linewidth=1.5)

    ax.set_xlabel(r"Photon index $\Gamma$", fontsize=12)
    ax.set_ylabel("Normalized distribution", fontsize=12)
    ax.set_title(f"Γ distribution by classification | {self.field_type}", fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.2, linestyle="--")

    plots_dir = self.base_dir / "Catalogues" / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    outpath = plots_dir / f"gamma_distribution_{self.field_type}.png"
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    self.logger.info(f"Γ distribution plot saved to {outpath}")




def post_analysis_plots(self):
    _plot_gamma_daox(self)
    _plot_gamma_counts(self)
    _plot_targets_on_lxluv(self)
    _plot_hr_gamma(self)
    _plot_gamma_distribution(self)





def _plot_gamma_err_vs_counts(self):
 
    csv_path = self.base_dir / "Datasets" / f"xray_weakness_comparison_{self.field_type}.csv"
    if not csv_path.exists():
        self.logger.info("xray_weakness_comparison csv not found, skipping Γ error vs counts plot")
        return
 
    gamma_types = ["pn", "combined", "mos"]
    markers = {"pn": "x", "combined": "s", "mos": "^"}
    colors = {"weak": "blue", "normal": "black", "strong": "red"}
 
    data = {(cls, t): {"dof": [], "err": []} for cls in colors for t in gamma_types}
 
    with open(csv_path, newline="") as f:
        for r in csv.DictReader(f):
            if r["classification"] not in colors:
                continue
            for t in gamma_types:
                gamma_col, gamma_lo_col, gamma_hi_col, dof_col = f"gamma_{t}", f"gamma_{t}_lo", f"gamma_{t}_hi", f"dof_{t}"
                if not r[gamma_col] or not r[dof_col]:
                    continue
                gamma = float(r[gamma_col])
                gamma_lo_val = float(r[gamma_lo_col])
                gamma_hi_val = float(r[gamma_hi_col])
                if gamma_lo_val == 0 and gamma_hi_val == 0:
                    continue
                if gamma_lo_val > gamma or gamma_hi_val < gamma:
                    continue
                dof = float(r[dof_col])
                g = data[(r["classification"], t)]
                g["dof"].append(dof)
                g["err"].append(gamma_hi_val - gamma_lo_val)
 
    fig, ax = plt.subplots(figsize=(8, 6))
    for cls in colors:
        for t in gamma_types:
            g = data[(cls, t)]
            if not g["dof"]:
                continue
            ax.scatter(g["dof"], g["err"], color=colors[cls], marker=markers[t], s=40)
 
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Fit degrees of freedom (proxy for counts)", fontsize=12)
    ax.set_ylabel(r"Γ error width (hi $-$ lo)", fontsize=12)
    ax.set_title(f"Γ uncertainty vs counts | {self.field_type}", fontsize=12)
 
    class_handles = [Line2D([0], [0], marker="o", color=colors[cls], linestyle="None", markersize=6, label=cls) for cls in colors]
    type_handles = [Line2D([0], [0], marker=markers[t], color="gray", linestyle="None", markersize=6, label=t) for t in gamma_types]
    class_legend = ax.legend(handles=class_handles, fontsize=10, loc="upper right")
    ax.add_artist(class_legend)
    ax.legend(handles=type_handles, fontsize=10, loc="lower left")
    ax.grid(True, alpha=0.2, linestyle="--")
 
    plots_dir = self.base_dir / "Catalogues" / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    outpath = plots_dir / f"gamma_err_vs_counts_{self.field_type}.png"
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    self.logger.info(f"Γ error vs counts plot saved to {outpath}")
 
 
def _plot_counts_vs_daox(self):
 
    csv_path = self.base_dir / "Datasets" / f"xray_weakness_comparison_{self.field_type}.csv"
    if not csv_path.exists():
        self.logger.info("xray_weakness_comparison csv not found, skipping counts vs Δαox plot")
        return
 
    colors = {"weak": "blue", "normal": "black", "strong": "red"}
    data = {cls: {"dof": [], "daox": []} for cls in colors}
 
    with open(csv_path, newline="") as f:
        for r in csv.DictReader(f):
            if r["classification"] not in colors:
                continue
            if not r["dof_combined"] or not r["delta_aox"]:
                continue
            data[r["classification"]]["dof"].append(float(r["dof_combined"]))
            data[r["classification"]]["daox"].append(float(r["delta_aox"]))
 
    fig, ax = plt.subplots(figsize=(8, 6))
    for cls in colors:
        g = data[cls]
        if not g["dof"]:
            continue
        ax.scatter(g["dof"], g["daox"], color=colors[cls], s=40, label=cls)
 
    ax.axhline(self.xray_weak_threshold, color="blue", linewidth=0.8, linestyle="--")
    ax.axhline(self.xray_strong_threshold, color="red", linewidth=0.8, linestyle="--")
    ax.set_xscale("log")
    ax.set_xlabel("Combined fit degrees of freedom (proxy for counts)", fontsize=12)
    ax.set_ylabel(r"$\Delta\alpha_{ox}$", fontsize=12)
    ax.set_title(f"Δαox vs counts | {self.field_type}", fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.2, linestyle="--")
 
    plots_dir = self.base_dir / "Catalogues" / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    outpath = plots_dir / f"counts_vs_daox_{self.field_type}.png"
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    self.logger.info(f"Counts vs Δαox plot saved to {outpath}")
 
 
def _plot_counts_vs_fvar(self):
 
    csv_path = self.base_dir / "Datasets" / f"xray_weakness_comparison_{self.field_type}.csv"
    if not csv_path.exists():
        self.logger.info("xray_weakness_comparison csv not found, skipping counts vs Fvar plot")
        return
 
    colors = {"weak": "blue", "normal": "black", "strong": "red"}
    measured = {cls: {"dof": [], "fvar": [], "err": []} for cls in colors}
    limits = {cls: {"dof": [], "fvar": []} for cls in colors}
 
    with open(csv_path, newline="") as f:
        for r in csv.DictReader(f):
            if r["classification"] not in colors:
                continue
            if not r["dof_combined"] or r["fvar"] == "" or r["fvar_err"] == "":
                continue
            dof = float(r["dof_combined"])
            fvar = float(r["fvar"])
            fvar_err = float(r["fvar_err"])
            if fvar == 0.0:
                limits[r["classification"]]["dof"].append(dof)
                limits[r["classification"]]["fvar"].append(fvar_err)
            else:
                measured[r["classification"]]["dof"].append(dof)
                measured[r["classification"]]["fvar"].append(fvar)
                measured[r["classification"]]["err"].append(fvar_err)
 
    fig, ax = plt.subplots(figsize=(8, 6))
    for cls in colors:
        m = measured[cls]
        if m["dof"]:
            ax.errorbar(m["dof"], m["fvar"], yerr=m["err"], fmt="o", color=colors[cls],
                        markersize=6, capsize=3)
        l = limits[cls]
        if l["dof"]:
            ax.scatter(l["dof"], l["fvar"], color=colors[cls], marker="v", s=30, alpha=0.5)
 
    ax.set_xscale("log")
    ax.set_xlabel("Combined fit degrees of freedom (proxy for counts)", fontsize=12)
    ax.set_ylabel(r"$F_{var}$", fontsize=12)
    ax.set_title(f"Fvar vs counts | {self.field_type}", fontsize=12)
 
    class_handles = [Line2D([0], [0], marker="o", color=colors[cls], linestyle="None", markersize=6, label=cls) for cls in colors]
    limit_handle = Line2D([0], [0], marker="v", color="gray", linestyle="None", markersize=6, alpha=0.5, label="upper limit")
    ax.legend(handles=class_handles + [limit_handle], fontsize=10)
    ax.grid(True, alpha=0.2, linestyle="--")
 
    plots_dir = self.base_dir / "Catalogues" / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    outpath = plots_dir / f"counts_vs_fvar_{self.field_type}.png"
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    self.logger.info(f"Counts vs Fvar plot saved to {outpath}")
 



def count_plots(self):
    _plot_gamma_err_vs_counts(self)
    _plot_counts_vs_daox(self)
    _plot_counts_vs_fvar(self)
