from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import json
import csv
import re
import sys
import math
import numpy as np
import pandas as pd

def _sdss_name_to_coords(name):
    stripped = name.replace("SDSS_J", "").replace("SDSS J", "")
    m = re.match(r"^(\d{2})(\d{2})(\d{2}\.?\d*)([+-]\d{2})(\d{2})(\d{2}\.?\d*)$", stripped)
    if not m:
        return None
    h, mnt, s, d, dm, ds = m.groups()
    ra = 15 * (int(h) + int(mnt) / 60 + float(s) / 3600)
    sign = -1 if d.startswith("-") else 1
    dec = sign * (abs(int(d)) + int(dm) / 60 + float(ds) / 3600)
    return ra, dec


def literature_check(self, match_radius_arcsec=5.0):
    fits_path = self.base_dir / "Catalogues" / "5XMM_DR15cat_source_v1.0.fits.gz"

    if not fits_path.exists():
        if self.auto_resolve:
            self.logger.info("5XMM catalogue not found, running Download_Catalogue")
            Download_Catalogue(self)
            if not fits_path.exists():
                self.logger.info("5XMM catalogue still missing, terminating")
                sys.exit(0)
        else:
            self.logger.info("5XMM catalogue not found, terminating")
            sys.exit(0)

    target_file = self.base_dir / "target_catalogue.json"
    with open(target_file) as f:
        targets = json.load(f)

    self.logger.info("Loading 5XMM-DR15 source catalogue...")
    with fits.open(fits_path) as hdul:
        cat = Table(hdul[1].data)
    self.logger.info(f"5XMM-DR15 catalogue loaded: {len(cat)} sources")

    cat_coords = SkyCoord(ra=cat["RA"] * u.deg, dec=cat["DEC"] * u.deg)

    results = []

    for name, data in targets.items():
        coords = _sdss_name_to_coords(name)
        if coords is not None:
            ra, dec = coords
        else:
            try:
                pos = SkyCoord.from_name(name)
                ra, dec = pos.ra.deg, pos.dec.deg
            except Exception as e:
                self.logger.info(f"{name}: could not resolve coordinates ({e}), skipping")
                continue

        target_coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
        sep = target_coord.separation(cat_coords).arcsec
        idx = sep.argmin()

        if sep[idx] > match_radius_arcsec:
            self.logger.info(f"{name}: no 5XMM match within {match_radius_arcsec}\" (closest: {sep[idx]:.1f}\")")
            continue

        row = {"target": name, "separation_arcsec": round(float(sep[idx]), 2)}
        for col in cat.colnames:
            row[col] = cat[idx][col]
        results.append(row)

        gamma = row.get("SPEC_GAMMA_PL", row.get("STACK_GAMMA", "N/A"))
        

    if results:
        out_path = self.base_dir / "Catalogues" / "5xmm_matched_targets.csv"
        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
            writer.writeheader()
            writer.writerows(results)
        self.logger.info(f"5XMM literature matches saved to {out_path}")
    else:
        self.logger.info("No 5XMM literature matches found")


def literature_consistency_check(self):
    lit_path = self.base_dir / "Catalogues" / "5xmm_matched_targets.csv"
    if not lit_path.exists():
        self.logger.info("5xmm_matched_targets.csv not found, run literature_check first")
        return

    comp_path = self.base_dir / "Datasets" / f"xray_weakness_comparison_{self.field_type}.csv"
    if not comp_path.exists():
        self.logger.info("xray_weakness_comparison csv not found, run xray_weakness_comparison first")
        return

    with open(lit_path, newline="") as f:
        lit_rows = {r["target"]: r for r in csv.DictReader(f)}

    def _best_gamma_for_target(rows_for_target, gamma_min, gamma_max):
        best = None
        for r in rows_for_target:
            for prefix in ("gamma_combined", "gamma_pn", "gamma_mos"):
                val = r.get(prefix)
                lo = r.get(f"{prefix}_lo")
                hi = r.get(f"{prefix}_hi")
                try:
                    val, lo, hi = float(val), float(lo), float(hi)
                except (TypeError, ValueError):
                    continue
                width = hi - lo
                pegged = val <= gamma_min or val >= gamma_max
                if width <= 0 or pegged:
                    continue
                if best is None or width < best[1]:
                    best = (prefix, width, val, lo, hi)
        return best

    rows_by_target = {}
    with open(comp_path, newline="") as f:
        for r in csv.DictReader(f):
            rows_by_target.setdefault(r["target"], []).append(r)

    own_rows = {}
    for target, rows_for_target in rows_by_target.items():
        best = _best_gamma_for_target(rows_for_target, self.specfit_gamma_min, self.specfit_gamma_max)
        if best is None:
            continue
        prefix, width, val, lo, hi = best
        own_rows[target] = {
            "gamma_combined": val if prefix == "gamma_combined" else "",
            "gamma_combined_lo": lo if prefix == "gamma_combined" else "",
            "gamma_combined_hi": hi if prefix == "gamma_combined" else "",
            "gamma_pn": val if prefix == "gamma_pn" else "",
            "gamma_pn_lo": lo if prefix == "gamma_pn" else "",
            "gamma_pn_hi": hi if prefix == "gamma_pn" else "",
            "gamma_mos": val if prefix == "gamma_mos" else "",
            "gamma_mos_lo": lo if prefix == "gamma_mos" else "",
            "gamma_mos_hi": hi if prefix == "gamma_mos" else "",
        }

    for target, lrow in lit_rows.items():
        orow = own_rows.get(target)
        if orow is None:
            continue

        sep = lrow.get("separation_arcsec", "N/A")

        g_own = orow.get("gamma_combined") or orow.get("gamma_pn") or orow.get("gamma_mos")
        g_own_lo = orow.get("gamma_combined_lo") or orow.get("gamma_pn_lo") or orow.get("gamma_mos_lo")
        g_own_hi = orow.get("gamma_combined_hi") or orow.get("gamma_pn_hi") or orow.get("gamma_mos_hi")
        g_own = float(g_own) if g_own not in ("", None) else None
        g_own_lo = float(g_own_lo) if g_own_lo not in ("", None) else 0.0
        g_own_hi = float(g_own_hi) if g_own_hi not in ("", None) else 0.0

        if g_own is None:
            continue

        g_lit = lrow.get("SPEC_GAMMA_PL")
        g_lit_lo = lrow.get("SPEC_GAMMA_ERR_LO_PL")
        g_lit_hi = lrow.get("SPEC_GAMMA_ERR_UP_PL")

        try:
            g_lit = float(g_lit)
            if math.isnan(g_lit):
                raise ValueError
        except (TypeError, ValueError):
            g_lit = lrow.get("STACK_GAMMA")
            try:
                g_lit = float(g_lit)
                if math.isnan(g_lit):
                    raise ValueError
                # STACK_GAMMA_ERR_LO/UP are +/- deltas, not absolute bounds (unlike
                # SPEC_GAMMA_ERR_LO/UP_PL) -- convert to absolute bounds here so the
                # (hi - lo) / 2 width calculation below is correct for either source.
                stack_err_lo = float(lrow.get("STACK_GAMMA_ERR_LO", 0) or 0)
                stack_err_hi = float(lrow.get("STACK_GAMMA_ERR_UP", 0) or 0)
                g_lit_lo = g_lit - stack_err_lo
                g_lit_hi = g_lit + stack_err_hi
            except (TypeError, ValueError):
                g_lit = None


                        
        g_own_lo = pd.to_numeric(g_own_lo, errors="coerce")
        g_own_hi = pd.to_numeric(g_own_hi, errors="coerce")

        g_own_err = (g_own_hi - g_own_lo) / 2



        g_lit_lo = pd.to_numeric(g_lit_lo, errors="coerce")
        g_lit_hi = pd.to_numeric(g_lit_hi, errors="coerce")

        g_lit_err = (g_lit_hi - g_lit_lo) / 2

        if g_lit is None:
            continue

        combined_err = math.sqrt(g_own_err**2 + g_lit_err**2)
        sigma = abs(g_own - g_lit) / combined_err if combined_err > 0 else float("inf")

        self.logger.info(
            f"{target}: matched with 5XMM, sep={sep}\", "
            f"Γ_own={g_own:.2f}, Γ_lit={g_lit:.2f}, consistent within {sigma:.2f}σ"
        )
