from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack
import astropy.units as u
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import sys

def _crossmatch(self):
    cat_dir = self.base_dir / "Catalogues"
    output_path = cat_dir / "matched_catalogue.fits"

    if output_path.exists():
        self.logger.info("Matched catalogue already exists, skipping crossmatch")
        return output_path

    self.logger.info("Loading 4XMM-DR14...")
    with fits.open(cat_dir / "4XMM_slim_DR14.fits.gz") as hdul:
        xmm = Table(hdul[1].data)

    self.logger.info("Loading SDSS DR16...")
    with fits.open(cat_dir / "SDSS_DR16_quasars.fits") as hdul:
        sdss = Table(hdul[1].data)

    self.logger.info(f"SDSS: {len(sdss)} quasars before redshift cut")
    sdss = sdss[(sdss["Z"] >= 1.8) & (sdss["Z"] <= 2.2)] 
    self.logger.info(f"SDSS: {len(sdss)} quasars after z = 1.8-2.2 cut")

    self.logger.info("Building sky coordinates...")
    xmm_coords = SkyCoord(ra=xmm["SC_RA"] * u.deg, dec=xmm["SC_DEC"] * u.deg)
    sdss_coords = SkyCoord(ra=sdss["RA"] * u.deg, dec=sdss["DEC"] * u.deg)

    self.logger.info("Crossmatching...")
    idx, sep, _ = sdss_coords.match_to_catalog_sky(xmm_coords)
    match_mask = sep < 5 * u.arcsec

    self.logger.info(f"Found {match_mask.sum()} matches within 5 arcsec")

    matched_sdss = sdss[match_mask]
    matched_xmm = xmm[idx[match_mask]]
    matched_sep = sep[match_mask]

    matched = hstack([matched_sdss, matched_xmm])
    matched["SEP_ARCSEC"] = matched_sep.to(u.arcsec).value

    matched.write(output_path, format="fits", overwrite=True)
    self.logger.info(f"Matched catalogue saved to {output_path.name}")

    return output_path




def _calculate_del_aox(self):
    from astropy.cosmology import FlatLambdaCDM

    cat_dir = self.base_dir / "Catalogues"
    input_path = cat_dir / "matched_catalogue.fits"
    output_path = cat_dir / f"aox_catalogue_{self.field_type}.fits"

    if output_path.exists():
        self.logger.info("Δαox catalogue already exists, skipping")
        return output_path

    if not input_path.exists():
        if self.auto_resolve:
            self.logger.info("Matched catalogue not found, running _crossmatch")
            _crossmatch(self)
            if not input_path.exists():
                self.logger.info("Matched catalogue still missing, terminating")
                sys.exit(0)
        else:
            self.logger.info("Matched catalogue not found, terminating")
            sys.exit(0)

    with fits.open(input_path) as hdul:
        t = Table(hdul[1].data)

    t = t[t["SC_SUM_FLAG"] < 3]
    t = t[t["SC_EP_3_FLUX"] > 0]
    t = t[t["SC_EP_4_FLUX"] > 0]
    t = t[t["PSFFLUX"][:, 3] > 0.1]
    t = t[t["BI_CIV"] == 0]
    t = t[(t["FIRST_MATCHED"] == 0) | (t["FIRST_FLUX"] < 1.0)]
    gi = t["PSFMAG"][:, 1] - t["PSFMAG"][:, 3]
    gi_median = np.median(gi)
    t = t[np.abs(gi - gi_median) < 0.1]

    self.logger.info(f"After quality cuts: {len(t)} sources")

    cosmo = FlatLambdaCDM(H0=self.H0, Om0=self.omega_m)
    z = np.asarray(t["Z"], dtype=float)
    d_L = cosmo.luminosity_distance(z).to(u.cm).value


    f3 = np.asarray(t["SC_EP_3_FLUX"], dtype=float)
    f4 = np.asarray(t["SC_EP_4_FLUX"], dtype=float)

    log_f3 = np.log10(np.abs(f3) + 1e-40)
    log_f4 = np.log10(np.abs(f4) + 1e-40)
    log_e3 = np.log10(1.5)
    log_e4 = np.log10(3.25)
    log_e2 = np.log10(2.0)

    log_f2kev = log_f3 + (log_f4 - log_f3) * (log_e2 - log_e3) / (log_e4 - log_e3)
    f2kev = (10 ** log_f2kev) / 2.418e17

    f_iband = np.asarray(t["PSFFLUX"][:, 3], dtype=float)
    f2500 = f_iband * 3.631e-29

    L2kev = 4 * np.pi * d_L**2 * f2kev * (1 + z) ** (self.photon_index - 2)
    L2500 = 4 * np.pi * d_L**2 * f2500 * (1 + z) ** (0.5 - 1)

    aox = 0.3838 * np.log10(L2kev / L2500)

    log_L2kev = np.log10(np.abs(L2kev) + 1e-40)
    log_L2500 = np.log10(np.abs(L2500) + 1e-40)

    valid = np.isfinite(log_L2kev) & np.isfinite(log_L2500) & np.isfinite(aox) & (f2kev > 0) & (f2500 > 0)

    coeffs = np.polyfit(log_L2500[valid], log_L2kev[valid], 1)
    log_L2kev_predicted = np.polyval(coeffs, log_L2500)

    delta_aox = log_L2kev - log_L2kev_predicted

    t["L2KEV"] = L2kev
    t["L2500"] = L2500
    t["AOX"] = aox
    t["DELTA_AOX"] = delta_aox
    t["VALID_AOX"] = valid.astype(int)

    t.write(output_path, format="fits", overwrite=True)
    self.logger.info(f"Δαox catalogue saved: {output_path.name}")

    return output_path



def Population_Analysis(self):
    _crossmatch(self)
    _calculate_del_aox(self)





def lxluv_plot(self, photon_index):

    cat_dir = self.base_dir / "Catalogues"
    input_path = cat_dir / f"aox_catalogue_{self.field_type}.fits"

    if not input_path.exists():
        if self.auto_resolve:
            self.logger.info("AOX catalogue not found, running _calculate_del_aox")
            _calculate_del_aox(self)
            if not input_path.exists():
                self.logger.info("AOX catalogue still missing, terminating")
                sys.exit(0)
        else:
            self.logger.info("AOX catalogue not found, terminating")
            sys.exit(0)

    with fits.open(input_path) as hdul:
        t = Table(hdul[1].data)

    valid = t["VALID_AOX"].data.astype(bool)
    tv = t[valid]

    log_L2500 = np.log10(tv["L2500"].data)
    log_L2kev = np.log10(tv["L2KEV"].data)
    delta_aox = tv["DELTA_AOX"].data

    xweak = delta_aox < -0.3

    coeffs = np.polyfit(log_L2500, log_L2kev, 1)
    x_fit = np.linspace(log_L2500.min(), log_L2500.max(), 200)
    y_fit = np.polyval(coeffs, x_fit)

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(log_L2500[~xweak], log_L2kev[~xweak],
               s=3, color="black", alpha=0.4, label="Normal")
    ax.scatter(log_L2500[xweak], log_L2kev[xweak],
               s=6, color="red", alpha=0.6, label="X-ray weak (Δαox < -0.3)")

    ax.plot(x_fit, y_fit, color="black", linewidth=1.0)
    ax.plot(x_fit, y_fit - 0.3, color="black", linewidth=0.8, linestyle="--")

    ax.set_xlabel(r"$\log(L_{2500\AA})$ [erg s$^{-1}$ Hz$^{-1}$]", fontsize=12)
    ax.set_ylabel(r"$\log(L_{2\,\mathrm{keV}})$ [erg s$^{-1}$ Hz$^{-1}$]", fontsize=12)
    ax.set_title(f"{self.target_name}  | {self.field_type} |  z = 1.8-2.2, Γ = {photon_index}", fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.2, linestyle="--")

    plots_dir = cat_dir / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    outpath = plots_dir / f"lx_luv_{self.field_type}_field_{photon_index}.png"
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    self.logger.info(f"LX-LUV plot saved to {outpath}")



def _civ_correlation(self, photon_index, civ_bins):
    cat_dir = self.base_dir / "Catalogues"
    input_path = cat_dir / f"aox_catalogue_{self.field_type}.fits"

    if not input_path.exists():
        if self.auto_resolve:
            self.logger.info("AOX catalogue not found, running _calculate_del_aox")
            _calculate_del_aox(self)
            if not input_path.exists():
                self.logger.info("AOX catalogue still missing, terminating")
                sys.exit(0)
        else:
            self.logger.info("AOX catalogue not found, terminating")
            sys.exit(0)

    with fits.open(input_path) as hdul:
        t = Table(hdul[1].data)

    valid = t["VALID_AOX"].data.astype(bool)
    tv = t[valid]

    c_kms = 2.998e5
    z = tv["Z"].data.astype(float)
    z_civ = tv["Z_CIV"].data.astype(float)
    delta_aox = tv["DELTA_AOX"].data.astype(float)

    blueshift = c_kms * (z - z_civ) / (1 + z)

    civ_valid = np.isfinite(blueshift) & np.isfinite(delta_aox) & (tv["ZWARN_CIV"].data == 0)
    bs = blueshift[civ_valid]
    da = delta_aox[civ_valid]
    physical = (bs > -1000) & (bs < 10000)
    bs = bs[physical]
    da = da[physical]

    
    bin_edges = np.percentile(bs, np.linspace(0, 100, civ_bins + 1))
    bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_medians = np.zeros(civ_bins)
    bin_errors = np.zeros(civ_bins)

    for i in range(civ_bins):
        mask = (bs >= bin_edges[i]) & (bs < bin_edges[i + 1])
        if mask.sum() > 0:
            bin_medians[i] = np.median(da[mask])
            bin_errors[i] = np.std(da[mask]) / np.sqrt(mask.sum())

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    axes[0].scatter(da, bs, s=2, color="black", alpha=0.3)
    axes[0].axhline(-0.3, color="red", linewidth=0.8, linestyle="--")
    axes[0].axhline(0, color="black", linewidth=0.6, linestyle="--")
    axes[0].set_ylabel("CIV blueshift (km/s)", fontsize=12)
    axes[0].set_xlabel(r"$\Delta\alpha_{ox}$", fontsize=12)
    axes[0].set_title("Δαox vs CIV blueshift ", fontsize=12)
    axes[0].grid(True, alpha=0.2, linestyle="--")

    axes[1].errorbar(bin_centres, bin_medians, yerr=bin_errors,
                     fmt="o", color="black", markersize=4,
                     elinewidth=0.8, capsize=3, capthick=0.8)
    axes[1].axhline(-0.3, color="red", linewidth=0.8, linestyle="--")
    axes[1].axhline(0, color="black", linewidth=0.6, linestyle="--")
    axes[1].set_xlabel("CIV blueshift (km/s)", fontsize=12)
    axes[1].set_ylabel(r"Median $\Delta\alpha_{ox}$", fontsize=12)
    axes[1].set_title("Binned median Δαox vs CIV blueshift", fontsize=12)
    axes[1].grid(True, alpha=0.2, linestyle="--")

    plt.suptitle(f"4XMM-DR14 × SDSS DR16  |  z = 1.8–2.2, Γ = {photon_index}  |  {self.field_type}", fontsize=12)
    plt.tight_layout()

    plots_dir = cat_dir / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    outpath = plots_dir / f"civ_correlation_{self.field_type}_{photon_index}.png"
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    self.logger.info(f"CIV correlation plot saved to {outpath}")



def _hubble_plot(self, photon_index):
    cat_dir = self.base_dir / "Catalogues"
    input_path = cat_dir / f"aox_catalogue_{self.field_type}.fits"

    if not input_path.exists():
        if self.auto_resolve:
            self.logger.info("AOX catalogue not found, running _calculate_del_aox")
            _calculate_del_aox(self)
            if not input_path.exists():
                self.logger.info("AOX catalogue still missing, terminating")
                sys.exit(0)
        else:
            self.logger.info("AOX catalogue not found, terminating")
            sys.exit(0)

    with fits.open(input_path) as hdul:
        t = Table(hdul[1].data)

    valid = t["VALID_AOX"].data.astype(bool)
    tv = t[valid]

    cosmo = FlatLambdaCDM(H0=self.H0, Om0=self.omega_m)

    z = tv["Z"].data.astype(float)
    f_x = tv["SC_EP_3_FLUX"].data.astype(float)
    f_uv = tv["PSFFLUX"][:, 3].data.astype(float) * 3.631e-29

    log_fx = np.log10(f_x)
    log_fuv = np.log10(f_uv)

    coeffs = np.polyfit(log_fuv, log_fx, 1)
    
    beta = coeffs[1]
    gamma = coeffs[0]
    log_dL_qso = (log_fx - photon_index * log_fuv - beta) / (2 - 2 * gamma)
    dL_qso = 10 ** log_dL_qso
    mu_qso = 5 * np.log10(dL_qso / 3.085e24) + 25

    z_grid = np.linspace(1.75, 2.25, 500)
    dL_lcdm = cosmo.luminosity_distance(z_grid).to(u.pc).value
    mu_lcdm = 5 * np.log10(dL_lcdm / 10)

    offset = np.median(mu_qso - 5 * np.log10(cosmo.luminosity_distance(z).to(u.pc).value / 10))

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.scatter(z, mu_qso - offset, s=2, color="black", alpha=0.3, label="Quasars")
    ax.plot(z_grid, mu_lcdm, color="red", linewidth=1.2, label=f"Flat ΛCDM (H0={self.H0}, Ωm={self.omega_m})")

    ax.set_xlabel("Redshift z", fontsize=12)
    ax.set_ylabel("Distance modulus μ", fontsize=12)
    ax.set_title(f"4XMM-DR14 × SDSS DR16  |  z = 1.8–2.2, Γ = {photon_index}  |  {self.field_type}", fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.2, linestyle="--")

    plots_dir = cat_dir / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    outpath = plots_dir / f"hubble_diagram_{self.field_type}_{photon_index}.png"
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    self.logger.info(f"Hubble diagram saved to {outpath}")


def civ_and_hubble(self, photon_index, civ_bins):
    _civ_correlation(self, photon_index, civ_bins)
    _hubble_plot(self, photon_index)
