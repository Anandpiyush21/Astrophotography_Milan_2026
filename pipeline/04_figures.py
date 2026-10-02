"""Report figures, the numbers quoted in the report, and the annotated image."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from PIL import Image
from scipy import ndimage as ndi
import config as C
from common import load_bands, load_fits, asinh_stretch, band_name

FIG = C.OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
FINAL = C.OUT / "final"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10})

prep = json.loads((C.WORK / "prep_stats.json").read_text())
comp = json.loads((C.WORK / "compose_stats.json").read_text())
cx, cy = prep["galaxy_centre_px"]
good = np.load(C.WORK / "coverage.npy")
skypix = np.load(C.WORK / "skypix.npy")
cube, _ = load_bands()
total = cube.sum(0)
exposure = load_fits("rel_exposure.fits")
L = load_fits("lum_adaptive.fits")
kernel = load_fits("scale_lum.fits")
sky = comp["sky_sigma"]


def show(ax, img, title, cmap="gray", **kw):
    ax.imshow(img, origin="lower", cmap=cmap, interpolation="nearest", **kw)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])


def cut(a, x0, x1, y0, y1):
    return a[y0:y1, x0:x1]


def save(name):
    plt.tight_layout(pad=1.5)
    plt.savefig(FIG / name, dpi=100)
    plt.close()


def bin16(a):
    return a.reshape(192, 16, 192, 16).sum((1, 3)).astype(float)


# fig 1: the raw bands
fig, axs = plt.subplots(2, 4, figsize=(14, 7.2))
h = 700
for i, ax in enumerate(axs.flat[:7]):
    d = cut(cube[i], cx - h, cx + h, cy - h, cy + h)
    d = d.reshape(h, 2, h, 2).sum((1, 3))
    show(ax, np.arcsinh(d), f"{band_name(C.BANDS[i])}\n{int(cube[i].sum()):,} counts", cmap="magma")
ax = axs.flat[7]
ax.bar(range(7), [c.sum() / 1e3 for c in cube], color=plt.cm.turbo(np.linspace(0.95, 0.1, 7)))
ax.set_xticks(range(7))
ax.set_xticklabels([f"{lo/1e3:g}-{hi/1e3:g}" for lo, hi in C.BANDS], rotation=45, fontsize=8)
ax.set_xlabel("band (keV)")
ax.set_ylabel("photons (thousands)")
ax.set_title("Photons per band")
save("fig1_raw_bands.png")

# fig 2: exposure correction
fig, axs = plt.subplots(1, 4, figsize=(16, 4.4))
tb, hb = bin16(total), bin16(cube[C.EXPOSURE_BAND])
show(axs[0], np.arcsinh(tb), "(a) all photons, 16x16 binned")
show(axs[1], ndi.gaussian_filter(hb, 0.7), "(b) 2.8-6 keV background", vmax=np.percentile(hb, 97))
im = axs[2].imshow(exposure, origin="lower", vmin=0, vmax=np.percentile(exposure[good], 99.5))
axs[2].set_title("(c) estimated relative exposure")
axs[2].set_xticks([])
axs[2].set_yticks([])
plt.colorbar(im, ax=axs[2], fraction=0.046)
eb = exposure[8::16, 8::16]
flat = np.where(eb > C.MIN_REL_EXPOSURE, tb / np.maximum(eb, 1e-3) - sum(prep["bkg_rates"]) * 256, 0)
show(axs[3], ndi.gaussian_filter(flat, 0.7), "(d) corrected, background removed", vmin=-10, vmax=60)
save("fig2_exposure.png")


def block_scatter(img, n=64):
    vals = [img[y:y + n, x:x + n][skypix[y:y + n, x:x + n]].mean()
            for y in range(0, 3072 - n, n) for x in range(0, 3072 - n, n)
            if skypix[y:y + n, x:x + n].mean() > 0.9]
    return 100 * float(np.std(vals) / abs(np.mean(vals)))


metrics = {"sky_scatter_raw": block_scatter(total),
           "sky_scatter_corrected": block_scatter(np.where(exposure > 0, total / np.maximum(exposure, 1e-3), 0))}

# fig 3: smoothing comparison on the NW wind
zoom = (cx - 60, cx + 540, cy + 100, cy + 700)
panels = [("lum_raw.fits", "(a) raw photons"),
          ("lum_gauss1.fits", "(b) Gaussian, 1 px"),
          ("lum_gauss8.fits", "(c) Gaussian, 8 px"),
          ("lum_adaptive.fits", "(d) adaptive")]
fig, axs = plt.subplots(1, 5, figsize=(20, 4.6))
for ax, (f, title) in zip(axs, panels):
    show(ax, asinh_stretch(cut(load_fits(f), *zoom), 1.5 * sky, 30), title, vmin=0, vmax=1)
im = axs[4].imshow(cut(kernel, *zoom), origin="lower", cmap="cividis_r", norm=LogNorm(0.5, 32))
axs[4].set_title("(e) kernel sigma used (px)")
axs[4].set_xticks([])
axs[4].set_yticks([])
plt.colorbar(im, ax=axs[4], fraction=0.046)
save("fig3_smoothing.png")


def fwhm(img, x, y, r=12):
    sub = img[y - r:y + r + 1, x - r:x + r + 1] - np.median(img[y - 3 * r:y + 3 * r, x - 3 * r:x + 3 * r])
    w = np.clip(sub, 0, None)
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    return float(2.3548 * np.sqrt((w * (xx ** 2 + yy ** 2)).sum() / w.sum() / 2))


yy, xx = np.mgrid[0:3072, 0:3072]
far = (np.hypot(yy - cy, xx - cx) > 400) & good
ps_y, ps_x = np.unravel_index(np.argmax(ndi.gaussian_filter(total, 1) * far), total.shape)
wind = np.zeros_like(good)
wind[cy + 250:cy + 450, cx + 50:cx + 250] = True
wind &= good & ~(ndi.gaussian_filter(total, 1.5) > 5 * ndi.gaussian_filter(total, 15))

metrics["point_source_xy"] = [int(ps_x), int(ps_y)]
metrics["smoothing"] = []
for f, name in [("lum_raw.fits", "raw"), ("lum_gauss1.fits", "Gaussian 1 px"),
                ("lum_gauss2.fits", "Gaussian 2 px"), ("lum_gauss4.fits", "Gaussian 4 px"),
                ("lum_gauss8.fits", "Gaussian 8 px"), ("lum_adaptive.fits", "adaptive")]:
    d = load_fits(f)
    rms = float(d[skypix].std())
    metrics["smoothing"].append({"method": name, "sky_rms": rms, "wind_snr": float(d[wind].mean() / rms),
                                 "fwhm": fwhm(d, ps_x, ps_y)})

# fig 4: the three colour groups
wide = comp["wide"]["box"]
names = {"R": "0.35-1.1 keV", "G": "1.1-2.2 keV", "B": "2.2-6 keV"}
fig, axs = plt.subplots(1, 4, figsize=(18, 5.4))
for ax, k, cmap in zip(axs, "RGB", ["inferno", "viridis", "cividis"]):
    d = cut(sum(load_fits(f"band{i}_colour.fits") for i in C.RGB_GROUPS[k]), *wide)
    show(ax, asinh_stretch(d, sky, np.percentile(d, 99.99)), names[k], cmap=cmap)
axs[3].imshow(Image.open(FINAL / "M82_Chandra_wide_Ramanujan.png"))
axs[3].set_title("combined")
axs[3].set_axis_off()
save("fig4_channels.png")

# fig 5: stretch and finishing on the core
core = comp["core"]
Lc = cut(L, *core["box"])
fig, axs = plt.subplots(1, 4, figsize=(18, 5.4))
show(axs[0], np.clip(Lc / core["white"], 0, 1), "(a) linear")
show(axs[1], asinh_stretch(Lc - C.BLACK_NSIG * sky, C.SOFT_NSIG["core"] * sky, core["white"]), "(b) asinh")
show(axs[2], np.load(C.WORK / "Ls_core.npy"), "(c) denoise, contrast, sharpen")
axs[3].imshow(Image.open(FINAL / "M82_Chandra_core_Ramanujan.png"))
axs[3].set_title("(d) with colour")
axs[3].set_axis_off()
save("fig5_postprocess.png")

# annotated version of the final image
img = np.asarray(Image.open(FINAL / "M82_Chandra_wide_Ramanujan.png"))
H, W = img.shape[:2]
x0, _, y0, _ = wide
dpi, fs = 200, 13
fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi)
ax = fig.add_axes([0, 0, 1, 1])
ax.imshow(img)
ax.set_axis_off()

bar = 2 * 60 / C.PIXEL_ARCSEC
kpc = 2 * C.DISTANCE_MPC * 1e3 * np.pi / (180 * 60)
ax.plot([W - 90 - bar, W - 90], [H - 90, H - 90], color="white", lw=3)
ax.text(W - 90 - bar / 2, H - 110, f"2′ ≈ {kpc:.1f} kpc", color="white", ha="center", fontsize=fs)

ox, oy, n = W - 140, 230, 110
arrow = dict(color="white", width=1.5, headwidth=8)
ax.annotate("", xy=(ox, oy - n), xytext=(ox, oy), arrowprops=arrow)
ax.annotate("", xy=(ox - n, oy), xytext=(ox, oy), arrowprops=arrow)
ax.text(ox, oy - n - 20, "N", color="white", ha="center", fontsize=fs)
ax.text(ox - n - 25, oy + 8, "E", color="white", ha="center", va="center", fontsize=fs)

for x, y, text, dx, dy in [(1543, 1525, "M82 X-1 / nucleus", 330, 60),
                           (cx + 120, cy + 450, "superwind", 280, -170),
                           (cx - 20, cy - 350, "southern outflow", -700, 250),
                           (1930, 2700, "the \"cap\"", -700, 30),
                           (cx - 75, cy + 30, "absorbed disk", -900, -40)]:
    X, Y = x - x0, H - (y - y0)
    ax.annotate(text, xy=(X, Y), xytext=(X + dx, Y + dy), color="white", fontsize=fs - 1,
                arrowprops=dict(arrowstyle="-", color="white", lw=0.8, alpha=0.8))

for j, k in enumerate("RGB"):
    ax.add_patch(Rectangle((60, 60 + j * 36), 26, 22, color=C.GROUP_COLOURS[k]))
    ax.text(100, 76 + j * 36, names[k], color="white", fontsize=fs - 1, va="center")
ax.text(60, H - 60, f"M82 · Chandra X-ray data: NASA/CXC · Processed by Hostel {C.HOSTEL}, IIT Hyderabad",
        color="white", fontsize=fs - 2, alpha=0.85)
plt.savefig(FINAL / "M82_Chandra_annotated_Ramanujan.png", dpi=dpi)
plt.close()

(C.WORK / "metrics.json").write_text(json.dumps(metrics, indent=2))
