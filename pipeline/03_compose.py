"""Combine the bands into the final colour image.

Brightness comes from all 7 bands summed. Colour comes from three energy groups
that were smoothed more heavily; only their ratios are used (Lupton et al. 2004).
"""
import json
import numpy as np
import tifffile
from PIL import Image
from scipy import ndimage as ndi
from skimage.restoration import denoise_nl_means
import config as C
from common import asinh_stretch, to_uint, save_png, load_fits

OUT = C.OUT / "final"
OUT.mkdir(parents=True, exist_ok=True)
cx, cy = json.loads((C.WORK / "prep_stats.json").read_text())["galaxy_centre_px"]
good = np.load(C.WORK / "coverage.npy")
sky = np.load(C.WORK / "skypix.npy")


def crop(a, box):
    x0, x1, y0, y1 = box
    return a[y0:y1, x0:x1]


def stretch(L, sky_sigma, white, soft, kernel):
    Ls = asinh_stretch(L - C.BLACK_NSIG * sky_sigma, soft * sky_sigma, white)

    # denoise the gas only; noise level measured on the gas, not the black sky
    hp = Ls - ndi.gaussian_filter(Ls, 2)
    sel = (Ls > 0.2) & (kernel > 1.5)
    noise = 1.4826 * np.median(np.abs(hp[sel] - np.median(hp[sel])))
    den = denoise_nl_means(Ls, patch_size=C.NLM["patch_size"], patch_distance=C.NLM["patch_distance"],
                           h=C.NLM["h_factor"] * noise, sigma=noise, fast_mode=True)
    detail = ndi.gaussian_filter(np.clip(C.DETAIL_KERNEL - kernel, 0, 1), 1)
    Ls = detail * Ls + (1 - detail) * den

    s, a = C.LOCAL_CONTRAST
    Ls += a * (Ls - ndi.gaussian_filter(Ls, s))
    s, a = C.SHARPEN
    Ls += a * detail * (Ls - ndi.gaussian_filter(Ls, s))
    return np.clip(Ls, 0, 1)


def colourise(Ls, chans, region):
    ch = np.stack([ndi.gaussian_filter(c, C.COLOUR_BLUR) for c in chans], -1)
    ch = (ch / ch[region].clip(0, None).sum(0)).clip(0, None)
    col = ch @ np.array([C.GROUP_COLOURS[k] for k in "RGB"])
    col /= col[region].sum(0) / col[region].sum()          # galaxy as a whole comes out neutral

    tot = col.sum(-1, keepdims=True)
    c = np.where(tot > 0, col / np.maximum(tot, 1e-12), 1 / 3)
    lo, hi = C.SAT_RAMP
    sat = 1 + (C.SATURATION - 1) * np.clip((Ls - lo) / (hi - lo), 0, 1)[..., None]
    c = (1 / 3 + sat * (c - 1 / 3)).clip(0, None)
    c /= c.sum(-1, keepdims=True)

    rgb = Ls[..., None] * 3 * c
    rgb /= np.maximum(rgb.max(-1, keepdims=True), 1)
    # push the brightest parts towards white instead of flat saturated colour
    h = np.clip((Ls - C.HIGHLIGHT_START) / (1 - C.HIGHLIGHT_START), 0, 1)[..., None] ** 1.5
    rgb += (Ls[..., None] - rgb) * h * (rgb.max(-1, keepdims=True) < 1.01)
    rgb += (1 - rgb) * h * 0.6
    return np.clip(rgb, 0, 1)


def export(rgb, stem):
    tifffile.imwrite(OUT / f"{stem}.tif", to_uint(rgb[::-1], 16), photometric="rgb", compression="zlib")
    save_png(OUT / f"{stem}.png", rgb)
    Image.fromarray(to_uint(rgb[::-1])).save(OUT / f"{stem}.jpg", quality=97, subsampling=0)


if __name__ == "__main__":
    L = load_fits("lum_adaptive.fits")
    kernel = load_fits("scale_lum.fits")
    bands = [load_fits(f"band{i}_colour.fits") for i in range(len(C.BANDS))]
    groups = [sum(bands[i] for i in C.RGB_GROUPS[k]) for k in "RGB"]
    sky_sigma = float(L[sky].std())
    yy, xx = np.mgrid[0:L.shape[0], 0:L.shape[1]]
    inner = np.hypot(yy - cy, xx - cx) < 600

    h = C.CORE_HALF
    boxes = {"wide": C.WIDE_BOX, "core": (cx - h, cx + h, cy - h, cy + h)}
    stats = {"sky_sigma": sky_sigma}
    for tag, box in boxes.items():
        Lc = crop(L, box)
        white = float(np.percentile(Lc[crop(good, box)], C.WHITE_PCT[tag]))
        Ls = stretch(Lc, sky_sigma, white, C.SOFT_NSIG[tag], crop(kernel, box))
        Ls *= np.clip(ndi.distance_transform_edt(crop(good, box)) / 40, 0, 1)   # fade out at mosaic edge
        export(colourise(Ls, [crop(g, box) for g in groups], crop(inner, box)), f"M82_Chandra_{tag}_Ramanujan")
        np.save(C.WORK / f"Ls_{tag}.npy", Ls)
        stats[tag] = {"box": [int(v) for v in box], "white": white}
    (C.WORK / "compose_stats.json").write_text(json.dumps(stats, indent=2))
