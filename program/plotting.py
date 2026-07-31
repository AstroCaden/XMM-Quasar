from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from program.processing import light_curve_processing, spectral_processing
import csv
import sys


def _plot_light_curves(self):
    for obsid in self.ObsID_current:
        fits_dir = self._ObsID_Paths(obsid)["fits_dir"]

        src_lc_files = list(fits_dir.glob("*_src_lc.ds"))
        if not src_lc_files:
            if self.auto_resolve:
                self.logger.info(f"{obsid}: source light curves not found, running light_curve_processing")
                light_curve_processing(self, self.energy_min_lc, self.energy_max_lc, self.bin_size_lc, self.CCF_path)
                src_lc_files = list(fits_dir.glob("*_src_lc.ds"))
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
            lc_file = fits_dir / f"{name}_src_lc.ds"
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
            lc_file = fits_dir / f"{name}_src_lc.ds"
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

            bin_size = 50
            n_bins = len(counts) // bin_size
            energy_binned = energy_kev[:n_bins * bin_size].reshape(n_bins, bin_size).mean(axis=1)
            counts_binned = counts[:n_bins * bin_size].reshape(n_bins, bin_size).sum(axis=1)
            errors_binned = np.sqrt(np.maximum(counts_binned, 1))

            mask = (counts_binned > 0) & (energy_binned >= 0.3) & (energy_binned <= 10.0)

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
            ax.set_xlim(0.3, 10.0)
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

        rows = []
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
                rows.append([name, e, c, err])

        with open(plots_dir / f"spectra_{obsid}.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["instrument", "energy_kev", "counts", "error"])
            writer.writerows(rows)

                

def plotting(self):
    _plot_light_curves(self)
    _plot_spectra(self)
