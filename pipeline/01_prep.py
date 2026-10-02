"""Coverage, relative exposure map, background rates and the data gap near X-1.

No exposure maps came with the data, so the exposure is estimated from the
2.8-6 keV band, where the particle background dominates away from the galaxy.
"""
import json
import warnings
import numpy as np
from scipy import ndimage as ndi
from astropy.io import fits
from skimage.registration import phase_cross_correlation
import config as C
from common import load_bands, save_fits, band_name

C.WORK.mkdir(exist_ok=True)
cube, hdr = load_bands()
ny, nx = cube.shape[1:]
total = cube.sum(0)
b = C.BKG_BIN


def binsum(a):
    return a.reshape(ny // b, b, nx // b, b).sum((1, 3))


smooth = ndi.gaussian_filter(total, 25)
cy, cx = np.unravel_index(np.argmax(smooth), smooth.shape)
yy, xx = np.mgrid[0:ny, 0:nx]
r = np.hypot(yy - cy, xx - cx)

cov_b = ndi.binary_fill_holes(ndi.binary_closing(binsum(total) > 0, iterations=2))
cov_b = ndi.binary_erosion(cov_b)
coverage = np.kron(cov_b, np.ones((b, b), bool))

# point sources: 5 sigma excess of a 1.5 px blur over a 20 px blur
n = 4 * np.pi * 1.5 ** 2
s1, s20 = ndi.gaussian_filter(total, 1.5), ndi.gaussian_filter(total, 20)
sig = (s1 - s20) * n / np.sqrt(np.maximum(s20 * n, 1))
sources = ndi.binary_dilation(sig > 5, iterations=6) | (r < C.GALAXY_MASK_R)

hard = binsum(cube[C.EXPOSURE_BAND]).astype(float)
masked = ndi.binary_dilation(binsum(sources.astype(int)) > 0)
hard[~cov_b | masked] = np.nan
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    bkg = ndi.generic_filter(hard, np.nanmedian, size=C.BKG_MEDIAN_BOX, mode="nearest")

# fill the hole under the galaxy with a Laplace fill (smooth, no streaks)
known = np.isfinite(bkg) & cov_b
hole = cov_b & ~known
bkg = np.where(known, bkg, np.nanmedian(bkg[known]))
for _ in range(4000):
    p = np.pad(bkg, 1, mode="edge")
    bkg[hole] = (0.25 * (p[:-2, 1:-1] + p[2:, 1:-1] + p[1:-1, :-2] + p[1:-1, 2:]))[hole]
bkg = ndi.median_filter(np.where(cov_b, bkg, 0), size=3)

exposure = ndi.zoom(bkg / np.median(bkg[cov_b & (bkg > 0)]), b, order=1).astype(np.float32)
good = coverage & (exposure >= C.MIN_REL_EXPOSURE)
exposure[~good] = 0

# zero-count strip inside bright emission next to M82 X-1
local = ndi.gaussian_filter(total, 2)
gap = (total == 0) & (local > 3) & good
gap = ndi.binary_dilation(gap) & (total < 0.2 * local)

# how unlikely the strip is as real sky: expected counts from the pixels around it
ok = (~gap).astype(float)
lam = ndi.gaussian_filter(total * ok, 4) / np.maximum(ndi.gaussian_filter(ok, 4), 1e-6)
lam_gap = lam[gap]

# the headers match, but check that the pixels really line up: each band against the
# all-band sum, central 1000 x 1000 px, 3 px pre-blur, 1/20 px precision
wcs_keys = ["CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CDELT1", "CDELT2"]
headers = [fits.getheader(C.RAW_DIR / f"m82_{lo}_{hi}_eV.fits") for lo, hi in C.BANDS]
wcs_same = all(all(h[k] == headers[0][k] for k in wcs_keys) for h in headers)
box = np.s_[cy - 500:cy + 500, cx - 500:cx + 500]
ref = ndi.gaussian_filter(total[box], 3)
shifts = [float(np.hypot(*phase_cross_correlation(ref, ndi.gaussian_filter(cube[i][box], 3),
                                                  upsample_factor=20)[0])) for i in range(len(C.BANDS))]

sky = good & ~sources & (r > C.BKG_FIT_R)
rates = [float(cube[i][sky].sum() / exposure[sky].sum()) for i in range(len(C.BANDS))]
inner = r < 600
py, px = np.unravel_index(np.argmax(cube[-1]), cube[-1].shape)

stats = {
    "galaxy_centre_px": [int(cx), int(cy)],
    "bkg_rates": rates,
    "gap_pixels": int(gap.sum()),
    "gap_expected_mean": float(lam_gap.mean()),
    "gap_expected_min": float(lam_gap.min()),
    "gap_expected_total": float(lam_gap.sum()),
    "wcs": {k: headers[0][k] for k in wcs_keys},
    "wcs_identical": bool(wcs_same),
    "registration_shift_px": shifts,
    "zero_frac_in_footprint": float((total[good] == 0).mean()),
    "total_photons": int(total.sum()),
    "brightest_pixel": {"x": int(px), "y": int(py), "counts": float(total[py, px]),
                        "hard_fraction": float(cube[-1][py, px] / total[py, px])},
    "bands": [{"band": band_name(band), "total_counts": int(cube[i].sum()),
               "net_counts_r600": float(cube[i][inner].sum() - rates[i] * exposure[inner].sum()),
               "bkg_rate_per_px": rates[i], "max_pixel": float(cube[i].max()),
               "nonzero_frac": float((cube[i] > 0).mean())}
              for i, band in enumerate(C.BANDS)],
}

np.save(C.WORK / "coverage.npy", good)
np.save(C.WORK / "skypix.npy", sky)
np.save(C.WORK / "gapmask.npy", gap)
save_fits("rel_exposure.fits", exposure, hdr)
(C.WORK / "prep_stats.json").write_text(json.dumps(stats, indent=2))
print(f"centre {cx},{cy}  gap pixels {gap.sum()}  background {sum(rates):.3f} cts/px")
