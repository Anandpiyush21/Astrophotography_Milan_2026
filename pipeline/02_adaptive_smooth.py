"""Adaptive smoothing, similar in spirit to CIAO csmooth.

Each pixel is smoothed with the smallest gaussian in C.SCALES whose counts are
at least tau sigma above the background. Flux is G*(D - B) / G*E, i.e. counts
and exposure are smoothed with the same kernel.
"""
import json
import numpy as np
from scipy import ndimage as ndi
import config as C
from common import load_bands, save_fits, load_fits, FFTSmoother, kernel_sum_sq

cube, hdr = load_bands()
exposure = load_fits("rel_exposure.fits").astype(float)
good = np.load(C.WORK / "coverage.npy")
gap = np.load(C.WORK / "gapmask.npy")
rates = np.array(json.loads((C.WORK / "prep_stats.json").read_text())["bkg_rates"])
K = len(C.SCALES)
var_factor = [kernel_sum_sq(s) for s in C.SCALES]

cube[:, ~good] = 0
weight = (good & ~gap).astype(float)
norm = ndi.gaussian_filter(weight, 2)
for band in cube:
    band[gap] = (ndi.gaussian_filter(band * weight, 2) / np.maximum(norm, 1e-6))[gap]

sm = FFTSmoother(exposure.shape)
F = sm.forward(exposure)
exp_k = [sm.smooth(F, s) for s in C.SCALES]


def scale_map(sig, tau):
    # fractional index of the first scale that reaches tau
    t = np.full(sig[0].shape, K - 1, np.float32)
    done = np.zeros(t.shape, bool)
    for k in range(K):
        hit = ~done & (sig[k] >= tau)
        if k == 0:
            t[hit] = 0
        else:
            frac = np.clip((tau - sig[k - 1]) / np.maximum(sig[k] - sig[k - 1], 1e-6), 0, 1)
            t[hit] = (k - 1 + frac)[hit]
        done |= hit
    # the median drops single noisy pixels that passed at a tiny scale
    t = ndi.median_filter(t, size=C.SCALEMAP_MEDIAN)
    return ndi.gaussian_filter(t, C.SCALEMAP_SMOOTH)


def smooth_image(counts, rate, taus):
    F = sm.forward(counts)
    d = [sm.smooth(F, s) for s in C.SCALES]
    sig = []
    for k in range(K):
        bk = rate * exp_k[k]
        sig.append((d[k] - bk) / np.sqrt(np.maximum(np.maximum(d[k], bk), 1e-9) * var_factor[k]))

    def flux(k):
        f = (d[k] - rate * exp_k[k]) / np.maximum(exp_k[k], 1e-3)
        f[~good] = 0
        return f
    return [scale_map(sig, tau) for tau in taus], flux


def blend(flux, t):
    out = np.zeros(t.shape, np.float32)
    for k in range(K):
        w = np.clip(1 - np.abs(t - k), 0, None)
        if w.any():
            out += w * flux(k)
    return out


total = cube.sum(0).astype(float)
(t_lum, t_col), flux = smooth_image(total, rates.sum(), [C.SNR_LUM, C.SNR_COLOUR])
save_fits("lum_adaptive.fits", blend(flux, t_lum), hdr)
save_fits("scale_lum.fits", np.interp(t_lum, np.arange(K), C.SCALES), hdr)

# fixed kernels, only for the comparison in the report
save_fits("lum_raw.fits", np.where(exposure > 0, (total - rates.sum() * exposure) / np.maximum(exposure, 1e-3), 0), hdr)
for s in (1.0, 2.0, 4.0, 8.0):
    save_fits(f"lum_gauss{s:g}.fits", flux(C.SCALES.index(s)), hdr)

for i, band in enumerate(cube):
    _, flux = smooth_image(band.astype(float), rates[i], [])
    save_fits(f"band{i}_colour.fits", blend(flux, t_col), hdr)
    print(f"band {i} done")
