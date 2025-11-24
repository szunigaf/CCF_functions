# CCF Functions

Tools to compute:
- Radial velocity (RV)
- Rotational velocity (v sin i)
- Bisector inverse slope (BIS)
- Bisector slope (b_b) and curvature (c_b)
- Anderson-Darling Gaussianity significance (ad_sig)

Supported Phase 3 instruments: ELODIE, UVES, UVES (combined arms), XSHOOTER, HARPS, NIRPS, FEROS.

The code was originally developed by P. Elliott and then extended by A. Bayo and S. Zúñiga-Fernández, and presented in [Zúñiga-Fernández et al. 2021](https://ui.adsabs.harvard.edu/abs/arXiv:2010.08575). It has evolved since publication for broader instrument support, configurability, and robustness with modern Python versions.

![CCF_functions](output_example-1.png)

---
## 1. Installation

Create or activate a Python 3.9+ environment and install dependencies:
```
pip install -r requirements.txt
```

Optional (interactive usage):
```
pip install jupyter ipywidgets
```

Directory expectations:
- `MODULES/` contains the core analysis script `calc_rv_new_chapman.py`
- `masks/` contains spectral masks (e.g. `K0.mas`, `M2_IR.mas`)
- `vsini_templates/` contains Gray rotational profile templates (`qqcc*.dat`)
- `INPUT/` (default) or user-specified directory contains FITS files
- `OUTPUTS/` is created automatically (stores CSV, PDF, and master summary)

---
## 2. Quick Start (Single File)

Example (HARPS spectrum in `HARPS/`):
```
python use_ccf_functions.py \
	--data_dir HARPS \
	--instrument HARPS \
	--mask K0 \
	--vel_min -120 --vel_max 120 \
	--vel_step 1.0 --wlen_step 0.02 --ncpu 8
```
Outputs:
- `OUTPUTS/master_output.txt` (appended summary line)
- `OUTPUTS/ccf_<object>_<mask>_<MJD>.csv` (raw CCF profile)
- `OUTPUTS/PDFS/ccf_<object>_<mask>_<MJD>.pdf` (diagnostic plot)

---
## 3. Batch Processing All Files in a Directory

Process all FEROS spectra in a folder:
```
python use_ccf_functions.py --data_dir FEROS --instrument FEROS --mask K0
```
Specify a pattern (default `*fits`):
```
python use_ccf_functions.py --data_dir HARPS_collection --file_pattern "ADP.*.fits" --instrument HARPS --mask R50G2
```

Resume mode skips already processed files (based on `master_output.txt`):
```
python use_ccf_functions.py --data_dir HARPS_collection --instrument HARPS --mask K0 --resume
```

---
## 4. UVES Combined Arms Mode

You can process two UVES arms as one spectrum.

Single pair:
```
python use_ccf_functions.py --instrument UVES_combine \
	--file1 UVES/ADP.2024-03-09T14:49:20.162.fits \
	--file2 UVES/ADP.2024-03-09T14:49:20.166.fits \
	--mask K0
```

Multiple pairs (pairs file with two paths per line):
```
python use_ccf_functions.py --instrument UVES_combine \
	--pairs_file uves_pairs.txt --mask K0 --resume
```
Example `uves_pairs.txt`:
```
# file1 file2
ADP.2024-03-09T14:49:20.162.fits ADP.2024-03-09T14:49:20.166.fits
...
```

Resume mode filters already processed pairs using a combined identifier `file1+file2` in `master_output.txt`.

---
## 5. Configuration File Usage

Instead of long CLI chains you can provide an INI file (see `ccf_config.ini`). Example:
```
[Basic]
data_dir = HARPS_collection
output_dir = OUTPUTS
mask = K0
instrument = HARPS

[CCF_Parameters]
vel_min = -120
vel_max = 120
vel_step = 1.0
wlen_step = 0.02
vsini_step = 1
ncpu = 8

[Processing]
file_pattern = ADP.*.fits
resume = true
verbose = true
```
Run:
```
python use_ccf_functions.py --config ccf_config.ini
```

Command-line arguments override config values if both supplied.

---
## 6. Jupyter Notebook Interactive Workflow

To do...

---
## 7. Parameter Guidelines

Velocity range (`vel_min`, `vel_max`):
- Typical: ±120–150 km/s unless high-velocity systems
- Narrower range speeds up processing.

`vel_step`:
- 1.0 km/s is a good balance
- <0.5 km/s sharply increases runtime.

`wlen_step`:
- 0.02 Å default; smaller captures finer structure at cost of speed
- <0.01 Å only if very high S/N and narrow lines.

`vsini_step`:
- 1 km/s standard; smaller improves precision for very slow rotators.

`ncpu`:
- Parallel chunking—diminishing returns >16 cores due to overhead.

`mask` selection:
- Match spectral type (e.g. K0 for solar-type, M2/M4 for late-type)
- IR masks (`M2_IR`, `M5_5_IR`) for NIRPS data if available.

---
## 8. Outputs Explained

`OUTPUTS/master_output.txt` columns:
```
file obj_name ra_deg dec_deg mjd_obs rv width depth vsini_calc_val min_chi snr_median mask central_wlen wlen_range BIS b_b b_b_err c_b ad_sig
```
- rv (km/s): Gaussian centroid
- width (sigma proxy): sqrt(width param)
- depth: Gaussian depth at centroid
- vsini_calc_val: best rotational broadening fit
- min_chi: residual measure from rotational fit
- BIS: bisector inverse slope (velocity span)
- b_b: bisector slope; b_b_err its std error
- c_b: curvature metric
- ad_sig: significance level from Anderson-Darling normality test

Individual CCF CSV: `ccf_<obj>_<mask>_<mjd>.csv` (velocity vs CCF profile).
Diagnostic plot: overlay of Gaussian fit, bisector panels, rotational profile comparison.

---
## 9. Resume Mode Logic

When `--resume` is used the script:
- Reads existing `master_output.txt`
- Extracts first token per line (file name or `file1+file2` for UVES_combine)
- Skips files/pairs already seen.

This prevents redundant re-processing during long runs.

---
## 10. Troubleshooting

Large CCF values:
- Often due to continuum normalization failure or mask mismatch.
- Check warnings printed: extreme flux range, problematic polynomial.

Mask not found:
- Ensure mask file exists in `masks/` and exact name matches choice.

Slow performance:
- Reduce velocity range or increase `vel_step`
- Increase `wlen_step` cautiously.

Empty or truncated output:
- Confirm FITS file has expected Phase 3 structure (WAVE, FLUX, ERR in HDU[1]).

UVES_combine errors:
- Verify both arms share object and have overlapping wavelength coverage.

`handle_normal_resume` not defined error (older versions):
- Update to latest version—function now included to filter already processed files.

LaTeX text errors (plot labels):
- If system lacks LaTeX, disable by removing `rc('text', usetex=True)` in `calc_rv_new_chapman.py`.

---
## 11. Citation

If you use this code in research, please cite:
```
Zúñiga-Fernández et al. 2021, A&A (original methodology)
```
And acknowledge the extended CCF toolkit (this repository).

---
## 12. Roadmap / Possible Extensions
- Extend the usage for other instruments
- Add jupyer-notebook usage example

Feel free to open issues or contribute improvements.

---
## 13. License

See `LICENSE` file for terms.

