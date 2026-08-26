from astropy.io import fits
import os
import subprocess
import shutil
from astropy.io import fits as astrofits

def _Instrume_Processing(self):
    unwanted_fits = []
    wanted_fits = []
    wanted_types = {"PN", "M1", "M2"}
    for obsid in self.ObsID_current:
        paths = self._ObsID_Paths(obsid)
        fits_dir = paths["fits_dir"]

        for file in fits_dir.glob("*.FIT"):
            name = file.name
            
            instrume = name[16:18]
            if instrume not in wanted_types:
                unwanted_fits.append(file)
            else:
                wanted_fits.append(file)
            
    return wanted_fits, unwanted_fits




                
def _obsid_to_fits(self):
    obsid_files = []
    for obsid in self.ObsID_current:
        paths = self._ObsID_Paths(obsid)
        obs_dir = paths["obs_dir"]
        fits_dir = paths["fits_dir"]
        for file in obs_dir.iterdir():
            if file.is_file():
                obsid_files.append(file.name)
                shutil.move(str(file), str(fits_dir / file.name))
    return obsid_files


def _set_ObsID(self, CCF_path):
    for obsid in self.ObsID_current:
        fits_dir = self._ObsID_Paths(obsid)["fits_dir"]
        os.environ["SAS_CCFPATH"] = CCF_path
        os.environ["SAS_CCF"] = str(fits_dir / "ccf.cif")
        for file in fits_dir.iterdir():
            if file.suffix == ".SAS" and "SUM" in file.name:
                os.environ["SAS_ODF"] = str(file)
                self.logger.info(f"SAS environment set: {file.name}")
             
       
def _run_SAS(self, CCF_path):
    for obsid in self.ObsID_current:
        paths = self._ObsID_Paths(obsid)
        fits_dir = paths["fits_dir"]

        existing_imaging = list(fits_dir.glob("*ImagingEvts.ds"))
        if existing_imaging:
            self.logger.info(f"{obsid}: SAS already run, skipping")
            continue

    

        for file in fits_dir.glob("*SUM.SAS"):
            file.unlink()
        for file in fits_dir.glob("ccf.cif"):
            file.unlink()

        os.environ["SAS_ODF"] = str(fits_dir)
        os.environ["SAS_CCFPATH"] = CCF_path

        subprocess.run(["cifbuild"], cwd=fits_dir, check=True)

        os.environ["SAS_CCF"] = str(fits_dir / "ccf.cif")

        subprocess.run(
            ["odfingest", f"odfdir={fits_dir}", "withodfdir=yes", f"outdir={fits_dir}"],
            cwd=fits_dir,
            check=True
        )

        sas_file = next(fits_dir.glob("*SUM.SAS"))
        os.environ["SAS_ODF"] = str(sas_file)
        
        pn_raw_files = list(fits_dir.glob("*PNU*ImagingEvts.ds")) + list(fits_dir.glob("*PNS*ImagingEvts.ds"))
        if pn_raw_files:
            pn_raw_file = pn_raw_files[0]
            subprocess.run(["epreject", f"eventset={pn_raw_file}"], cwd=fits_dir, check=True)
        subprocess.run(["epproc"], cwd=fits_dir, check=True)
        subprocess.run(["emproc"], cwd=fits_dir, check=True)








def _gti_filtering(self):
    for obsid in self.ObsID_current:
        fits_dir = self._ObsID_Paths(obsid)["fits_dir"]


        existing_clean = list(fits_dir.glob("*_clean.ds"))
        if existing_clean:
            self.logger.info(f"{obsid}: GTI filtering already done, skipping")
            continue

        imaging_files = list(fits_dir.glob("*ImagingEvts.ds"))
        
        if not imaging_files:
            if self.auto_resolve:
                self.logger.info(f"{obsid}: ImagingEvts not found, running _obsid_to_fits and _run_SAS")
                _obsid_to_fits(self)
                _run_SAS(self, self.CCF_path)
                imaging_files = list(fits_dir.glob("*ImagingEvts.ds"))
                if not imaging_files:
                    self.logger.info(f"{obsid}: ImagingEvts still missing after _run_SAS, skipping")
                    continue
            else:
                self.logger.info(f"{obsid}: ImagingEvts not found, terminating")
                sys.exit(0)
        

        pn_file = next(fits_dir.glob("*EPN*ImagingEvts.ds"), None)
        mos1_file = next(fits_dir.glob("*EMOS1*ImagingEvts.ds"), None)
        mos2_file = next(fits_dir.glob("*EMOS2*ImagingEvts.ds"), None)

        instruments = []
        if pn_file:
            instruments.append((pn_file, "pn", "10000:12000"))
        if mos1_file:
            instruments.append((mos1_file, "mos1", "10000:15000"))
        if mos2_file:
            instruments.append((mos2_file, "mos2", "10000:15000"))


        for file, name, energy in instruments:
            lc_file = fits_dir / f"{name}_bkg_lc.ds"
            gti_file = fits_dir / f"{name}_gti.ds"
            file_clean = fits_dir / f"{name}_clean.ds"

            subprocess.run([
                "evselect",
                f"table={file}",
                "withrateset=yes",
                f"rateset={lc_file}",
                "maketimecolumn=yes",
                "timebinsize=100",
                "makeratecolumn=yes",
                f"expression=#XMMEA_EP && (PI in [{energy}]) && (PATTERN==0)" if name == "pn"
                else f"expression=#XMMEA_EM && (PI in [{energy}])",
            ], cwd=fits_dir, check=True)

           
            threshold = "0.4" if name == "pn" else "0.35"
            subprocess.run([
                "tabgtigen",
                f"table={lc_file}",
                f"gtiset={gti_file}",
                f"expression=RATE<={threshold}",
            ], cwd=fits_dir, check=True)

           
            subprocess.run([
                "evselect",
                f"table={file}",
                "withfilteredset=yes",
                f"filteredset={file_clean}",
                f"expression=(gti({gti_file},TIME))",
                "keepfilteroutput=yes",
            ], cwd=fits_dir, check=True)

            self.logger.info("GTI filtering done")



def _sky_coordinates(self, fits_dir):
        for name in ["pn_clean.ds", "mos1_clean.ds", "mos2_clean.ds"]:
            candidate = fits_dir / name
            
            if candidate.exists():
                with fits.open(candidate) as hdul:
                    hdr = hdul["EVENTS"].header
                    ra_ref  = hdr["REFXCRVL"]
                    dec_ref = hdr["REFYCRVL"]
                    x_ref   = hdr["REFXCRPX"]
                    y_ref   = hdr["REFYCRPX"]
                    x_scale = hdr["REFXCDLT"]
                    y_scale = hdr["REFYCDLT"]
                    
                src_x = int(round(x_ref + (self.RA  - ra_ref) / x_scale))
                src_y = int(round(y_ref + (self.Dec - dec_ref) / y_scale))
                
                return src_x, src_y





        
def pre_processing(self, CCF_path):
    _Instrume_Processing(self)
    _set_ObsID(self, CCF_path)
    _obsid_to_fits(self)
    _run_SAS(self, CCF_path)
    _gti_filtering(self)





def spectral_processing(self, channel_min_spec, channel_max_spec, bin_size_spec, CCF_path):
    _set_ObsID(self, CCF_path)
    for obsid in self.ObsID_current:
        fits_dir = self._ObsID_Paths(obsid)["fits_dir"]

        if list(fits_dir.glob("*_src_spec.fits")):
            self.logger.info(f"{obsid}: spectral processing already done, skipping")
            continue

        clean_files = list(fits_dir.glob("*_clean.ds"))
        if not clean_files:
            if self.auto_resolve:
                self.logger.info(f"{obsid}: clean files not found, running _gti_filtering")
                _gti_filtering(self)
                clean_files = list(fits_dir.glob("*_clean.ds"))
                if not clean_files:
                    self.logger.info(f"{obsid}: clean files still missing after _gti_filtering, skipping")
                    continue
            else:
                self.logger.info(f"{obsid}: clean files not found, terminating")
                sys.exit(0)


        src_x, src_y = _sky_coordinates(self, fits_dir)

       
        src_r    = 600
        bkg_r_in = 1200
        bkg_r_out = 2400

        src_region = f"(X,Y) IN circle({src_x},{src_y},{src_r})"
        bkg_region = f"(X,Y) IN annulus({src_x},{src_y},{bkg_r_in},{bkg_r_out})"

        all_instruments = [
            ("pn",   fits_dir / "pn_clean.ds",   "#XMMEA_EP&&(PATTERN<=4)",  channel_min_spec, channel_max_spec),
            ("mos1", fits_dir / "mos1_clean.ds", "#XMMEA_EM&&(PATTERN<=12)", channel_min_spec, 11999),
            ("mos2", fits_dir / "mos2_clean.ds", "#XMMEA_EM&&(PATTERN<=12)", channel_min_spec, 11999),
        ]
        instruments = [(name, evt, flag, ch_min, ch_max)
                       for name, evt, flag, ch_min, ch_max in all_instruments
                       if evt.exists()]

        for f in fits_dir.glob("*_src_spec.fits"):
            f.unlink()
        for f in fits_dir.glob("*_bkg_spec.fits"):
            f.unlink()


        for name, evt, flag, ch_min, ch_max in instruments:
            src_spec = fits_dir / f"{name}_src_spec.fits"
            bkg_spec = fits_dir / f"{name}_bkg_spec.fits"
            rmf      = fits_dir / f"{name}.rmf"
            arf      = fits_dir / f"{name}.arf"

            subprocess.run([
                "evselect", f"table={evt}",
                "withspectrumset=yes", f"spectrumset={src_spec}",
                "energycolumn=PI",
                f"specchannelmin={ch_min}",
                f"specchannelmax={ch_max}",
                f"spectralbinsize={bin_size_spec}",
                "withspecranges=yes",
                "writedss=yes",
                f"expression={flag}&&{src_region}",
            ], cwd=fits_dir, check=True)

            subprocess.run([
                "evselect", f"table={evt}",
                "withspectrumset=yes", f"spectrumset={bkg_spec}",
                "energycolumn=PI",
                f"specchannelmin={ch_min}",
                f"specchannelmax={ch_max}",
                f"spectralbinsize={bin_size_spec}",
                "withspecranges=yes",
                "writedss=yes",
                f"expression={flag}&&{bkg_region}",
            ], cwd=fits_dir, check=True)

        
            subprocess.run(["backscale", f"spectrumset={src_spec}", f"badpixlocation={evt}"], cwd=fits_dir, check=True)
            subprocess.run(["backscale", f"spectrumset={bkg_spec}", f"badpixlocation={evt}"], cwd=fits_dir, check=True)

            
            subprocess.run(["rmfgen", f"spectrumset={src_spec}", f"rmfset={rmf}"], cwd=fits_dir, check=True)
            subprocess.run(["arfgen", f"spectrumset={src_spec}", f"arfset={arf}", f"withrmfset=yes", f"rmfset={rmf}", f"badpixlocation={evt}"], cwd=fits_dir, check=True)

            self.logger.info("Spectral processing done.")






def light_curve_processing(self, energy_min_lc, energy_max_lc, bin_size_lc, CCF_path):
    _set_ObsID(self, CCF_path)
    for obsid in self.ObsID_current:
        fits_dir = self._ObsID_Paths(obsid)["fits_dir"]


        if list(fits_dir.glob("*_corrected_lc.ds")):
            self.logger.info(f"{obsid}: light curve processing already done, skipping")
            continue

        clean_files = list(fits_dir.glob("*_clean.ds"))
        if not clean_files:
            if self.auto_resolve:
                self.logger.info(f"{obsid}: clean files not found, running _gti_filtering")
                _gti_filtering(self)
                clean_files = list(fits_dir.glob("*_clean.ds"))
                if not clean_files:
                    self.logger.info(f"{obsid}: clean files still missing after _gti_filtering, skipping")
                    continue
            else:
                self.logger.info(f"{obsid}: clean files not found, terminating")
                sys.exit(0)


        src_x, src_y = _sky_coordinates(self, fits_dir)

        src_r    = 600
        bkg_r_in = 1200
        bkg_r_out = 2400

        src_region = f"(X,Y) IN circle({src_x},{src_y},{src_r})"
        bkg_region = f"(X,Y) IN annulus({src_x},{src_y},{bkg_r_in},{bkg_r_out})"
        energy     = f"{energy_min_lc}:{energy_max_lc}"

        all_instruments = [
            ("pn",   fits_dir / "pn_clean.ds",   "#XMMEA_EP&&(PATTERN<=4)"),
            ("mos1", fits_dir / "mos1_clean.ds", "#XMMEA_EM&&(PATTERN<=12)"),
            ("mos2", fits_dir / "mos2_clean.ds", "#XMMEA_EM&&(PATTERN<=12)"),
        ]
        instruments = [(name, evt, flag)
                       for name, evt, flag in all_instruments
                       if evt.exists()]

        for name, evt, flag in instruments:
            src_lc = fits_dir / f"{name}_src_lc.ds"
            bkg_lc = fits_dir / f"{name}_bkg_lc_sci.ds"
            corrected_lc = fits_dir / f"{name}_corrected_lc.ds"

       
            subprocess.run([
                "evselect", f"table={evt}",
                "withrateset=yes", f"rateset={src_lc}",
                "maketimecolumn=yes", "makeratecolumn=yes",
                f"timebinsize={bin_size_lc}",
                f"expression={flag}&&(PI in [{energy}])&&{src_region}",
            ], cwd=fits_dir, check=True)

          
            subprocess.run([
                "evselect", f"table={evt}",
                "withrateset=yes", f"rateset={bkg_lc}",
                "maketimecolumn=yes", "makeratecolumn=yes",
                f"timebinsize={bin_size_lc}",
                f"expression={flag}&&(PI in [{energy}])&&{bkg_region}",
            ], cwd=fits_dir, check=True)



            try:
                subprocess.run([
                    "epiclccorr",
                    f"eventlist={evt}",
                    f"outset={corrected_lc}",
                    f"srctslist={src_lc}",
                    "applyabsolutecorrections=yes",
                    "withbkgset=yes",
                    f"bkgtslist={bkg_lc}",
                ], cwd=fits_dir, check=True)
            except subprocess.CalledProcessError:
                self.logger.info(f"{name}: epiclccorr failed, likely empty light curve, skipping")
                continue
