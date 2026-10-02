"""Builds the PDF report from the numbers the pipeline saved."""
import json
import datetime
import platform
import numpy as np
import astropy, scipy, skimage, matplotlib, PIL, weasyprint
from weasyprint import HTML
import config as C

prep = json.loads((C.WORK / "prep_stats.json").read_text())
met = json.loads((C.WORK / "metrics.json").read_text())
sm = {r["method"]: r for r in met["smoothing"]}
raw, g8, ad = sm["raw"], sm["Gaussian 8 px"], sm["adaptive"]
bp = prep["brightest_pixel"]
team = ", ".join(C.TEAM) if C.TEAM else "-"
F, FI = "../outputs/figures/", "../outputs/final/"
REPO = "github.com/Anandpiyush21/Astrophotography_Milan_2026"


def rows(items):
    return "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in items)


content = ["O VII / O VIII lines, coolest wind gas", "O and Fe-L lines, soft thermal wind",
           "Fe-L and Ne lines, brightest band", "Mg and Ne lines, warm gas", "Si lines, hot gas near the starburst",
           "S lines and hard continuum", "hard continuum: X-ray binaries, absorbed nucleus"]
band_table = rows((b["band"], f"{b['total_counts']:,}", f"{100 * b['total_counts'] / prep['total_photons']:.1f}%",
                   f"{b['net_counts_r600']:,.0f}",
                   f"{b['bkg_rate_per_px']:.4f}", int(b["max_pixel"]), f"{100 * b['nonzero_frac']:.2f}%", c)
                  for b, c in zip(prep["bands"], content))
wcs = prep["wcs"]
max_shift = max(prep["registration_shift_px"])
smooth_table = rows((r["method"], f"{r['sky_rms']:.3f}", f"{r['wind_snr']:.2f}", f"{r['fwhm']:.1f}")
                    for r in met["smoothing"])
skip = {"ROOT", "RAW_DIR", "WORK", "OUT", "REPORT", "TEAM", "HOSTEL"}
param_table = rows((f"<code>{k}</code>", f"<code>{getattr(C, k)}</code>")
                   for k in dir(C) if k.isupper() and k not in skip)
versions = rows([("Python", platform.python_version(), "pipeline language"),
                 ("NumPy", np.__version__, "array maths"),
                 ("SciPy", scipy.__version__, "Gaussian and median filters, FFTs, morphology"),
                 ("Astropy", astropy.__version__, "FITS input and output, header checks"),
                 ("scikit-image", skimage.__version__, "registration check, non-local-means denoising"),
                 ("Matplotlib", matplotlib.__version__, "figures in this report"),
                 ("Pillow", PIL.__version__, "PNG, JPEG and 16-bit TIFF export"),
                 ("WeasyPrint", weasyprint.__version__, "this PDF")])

html = f"""<!doctype html><html><head><meta charset="utf-8"><title>M82 report, Hostel {C.HOSTEL}</title>
<style>
@page {{ size: A4; margin: 16mm 15mm;
  @bottom-center {{ content: "Hostel {C.HOSTEL} · M82 · " counter(page) "/" counter(pages); font-size: 8pt; color: #888; }} }}
body {{ font-family: "DejaVu Sans", sans-serif; font-size: 9.6pt; line-height: 1.45; color: #222; }}
h1 {{ font-size: 20pt; margin: 0 0 1mm; }}
h2 {{ font-size: 13pt; border-bottom: 1px solid #bbb; margin-top: 7mm; padding-bottom: 1mm; }}
figure {{ margin: 3mm 0; page-break-inside: avoid; }}
figure img {{ width: 100%; }}
figcaption {{ font-size: 8.4pt; color: #555; }}
table.t {{ border-collapse: collapse; width: 100%; font-size: 8.6pt; margin: 2mm 0; page-break-inside: avoid; }}
table.t th {{ background: #333; color: #fff; text-align: left; padding: 1.2mm; }}
table.t td {{ border-bottom: 0.5px solid #ccc; padding: 1mm; }}
table.t td:first-child {{ white-space: nowrap; }}
.eq {{ text-align: center; font-family: "DejaVu Serif", serif; margin: 2mm 0; }}
code {{ font-size: 8pt; }}
.pb {{ page-break-before: always; }}
a {{ color: #1a5fb4; text-decoration: none; }}
pre.cmd {{ background: #f2f2f2; border-left: 3px solid #999; padding: 2mm 3mm; font-size: 8.4pt; line-height: 1.5; margin: 2mm 0; }}
</style></head><body>

<h1>M82 in X-rays</h1>
<p>Cepheid Astro-Photography Competition, IIT Hyderabad · stacking report<br>
Hostel: <b>{C.HOSTEL}</b> · Team: {team} · {datetime.date.today():%d %B %Y}<br>
Code: <a href="https://{REPO}">{REPO}</a></p>

<figure><img src="{FI}M82_Chandra_annotated_Ramanujan.png">
<figcaption>Final image, annotated. 18′ × 18′, north up, east left, 0.492″ per pixel.
Orange is 0.35–1.1 keV, green 1.1–2.2 keV, blue 2.2–6 keV. The submitted file is the clean
version, <code>M82_Chandra_wide_Ramanujan.png</code> (a 16-bit TIFF master is also kept). Data: NASA/CXC.</figcaption></figure>

<h2>Summary</h2>
<p>The data we were given are not a set of exposures to align and average. They are seven Chandra
X-ray images of M82, one per energy band, already on the same pixel grid. Most pixels hold zero or one
photon. For us, stacking meant two things: combining the seven bands into one image, and dealing with
the photon noise without blurring the point sources.</p>
<p>We did this in Python in five steps: estimate the missing exposure map, smooth each pixel only as much
as it needs, build brightness from all bands and colour from three energy groups, then stretch and finish.
In the faint wind the signal-to-noise per pixel went from {raw['wind_snr']:.2f} to {ad['wind_snr']:.1f}, while point sources
stayed as sharp as in the raw data ({ad['fwhm']:.1f} px against {raw['fwhm']:.1f} px). An ordinary 8 px blur gets
{g8['wind_snr']:.1f} and spreads the same source to {g8['fwhm']:.1f} px.</p>
<p>The whole process is scripted. A single command rebuilds the submitted image and this report from the raw
FITS files (Section 10). No manual painting, masking or retouching was used.</p>

<h2>1. M82</h2>
<p>M82 is a starburst galaxy about {C.DISTANCE_MPC} Mpc (12 million light-years) away, seen nearly edge-on.
A close pass by M81 set off very fast star formation in its centre. The supernovae that followed heat gas to
millions of degrees and push it out above and below the disk. This outflow, the superwind, is bright in soft
X-rays and is the orange structure in our image. The image also shows:</p>
<ul>
<li>point sources, mostly X-ray binaries, including the very bright M82 X-1 and, in the nucleus, the
ultraluminous pulsar M82 X-2;</li>
<li>the disk itself, which shows up blue because its gas and dust absorb the low-energy X-rays;</li>
<li>the "cap", a faint cloud about 11 kpc north of the galaxy where the wind runs into surrounding gas.</li>
</ul>

<h2>2. The data</h2>
<table class="t"><tr><th>Band</th><th>Photons</th><th>Share</th><th>Net photons (r &lt; 4.9′)</th><th>Background (cts/px)</th><th>Max pixel</th><th>Pixels with a photon</th><th>Main content</th></tr>{band_table}
<tr><td><b>Total</b></td><td><b>{prep['total_photons']:,}</b></td><td>100%</td><td></td><td></td><td></td><td></td><td></td></tr></table>
<ul>
<li>Each file is 3072 × 3072 pixels of 32-bit floats holding whole photon counts, at 0.492″ per pixel
(CDELT = {wcs['CDELT2']} deg), Chandra's native pixel size. The projection is tangent-plane (TAN), centred on
RA {wcs['CRVAL1']}°, Dec {wcs['CRVAL2']}°.</li>
<li>Even after summing all seven bands, {100 * prep['zero_frac_in_footprint']:.0f}% of the pixels inside the observed area hold no
photon at all. The noise is Poisson (shot) noise, not the Gaussian noise of a camera.</li>
<li>The field is a mosaic of several pointings taken at different angles, and no exposure maps were included.
Because of that, the detector chip edges show up as steps in the background (Fig. 2a).</li>
<li>The brightest pixel (x={bp['x']}, y={bp['y']}) has {bp['counts']:.0f} photons, {100 * bp['hard_fraction']:.0f}% of them
above 2.8 keV. This is M82 X-1, and it suffers from pile-up. Right next to it is a diagonal strip of
{prep['gap_pixels']} pixels with zero counts in every band. From the pixels around it, the strip should hold about
{prep['gap_expected_mean']:.0f} photons per pixel, {prep['gap_expected_total']:,.0f} in total. For Poisson counts the chance of getting none
is e<sup>−{prep['gap_expected_total']:,.0f}</sup>, so this is not real sky. It is most likely a column removed during the original processing.</li>
</ul>
<figure><img src="{F}fig1_raw_bands.png"><figcaption>Fig. 1. The seven raw bands (2×2 binned, asinh).
The wind dominates below 1.1 keV, and point sources and the nucleus dominate above 2.2 keV.</figcaption></figure>

<h2>3. Alignment and stacking</h2>
<p><b>Alignment check.</b> We compared the coordinate keys (CTYPE, CRVAL, CRPIX, CDELT) of all seven headers
key by key, and they are {"identical" if prep['wcs_identical'] else "<b>not</b> identical"}. Identical headers do not prove the pixels line up, so each band
was also cross-correlated against the all-band sum (phase correlation, central 1000 × 1000 px, 3 px pre-blur,
1/20 px precision). The largest offset was {max_shift:.2f} px, far below Chandra's 0.5″ resolution. We therefore did
not resample any band, which avoids interpolation blur.</p>
<p><b>Why a plain sum.</b> The pixels are photon counts. A sum of Poisson counts is again Poisson, and every
photon gets the same weight. A mean has the same S/N, and a median throws information away, because most
pixels are zero. So the brightness layer is the plain sum of the seven bands ({prep['total_photons']:,} photons), and
the bands are also kept apart for colour.</p>

<h2>4. Step 1: exposure map and background (<code>01_prep.py</code>)</h2>
<p>Counts are flux times exposure, so parts of the mosaic covered by more pointings look brighter.
Without a correction, these steps turn into visible seams once the image is smoothed.</p>
<p>The exposure had to be estimated from the data. In the 2.8–6 keV band, away from sources, most events
come from the detector's particle background, which builds up at a steady rate per second of exposure.
There the galaxy is also compact and easy to mask. So we:</p>
<ol>
<li>masked point sources and a {C.GALAXY_MASK_R} px circle around the nucleus;</li>
<li>binned the hard band to {C.BKG_BIN}×{C.BKG_BIN} px and took a median over {C.BKG_MEDIAN_BOX}×{C.BKG_MEDIAN_BOX} bins, which keeps the chip edges sharp;</li>
<li>filled the masked area with a Laplace (harmonic) fill and scaled the result to a median of 1 (Fig. 2c);</li>
<li>measured the background level of each band only beyond {C.BKG_FIT_R} px from the centre, so the wind is not
counted as background.</li>
</ol>
<p>The spread of the sky level between 64×64 px blocks fell from {met['sky_scatter_raw']:.0f}% to
{met['sky_scatter_corrected']:.0f}%, and the seams are gone (Fig. 2d). Pixels in the zero-count strip near X-1
were filled in from their neighbours.</p>
<figure><img src="{F}fig2_exposure.png"><figcaption>Fig. 2. (a) All photons, showing the mosaic.
(b) The hard-band background. (c) Estimated exposure. (d) After correction.</figcaption></figure>

<h2>5. Step 2: adaptive smoothing (<code>02_adaptive_smooth.py</code>)</h2>
<p>No single blur works for both the point sources and the faint wind. We used the approach of CIAO's
<code>csmooth</code> (Ebeling et al. 2006). For each pixel we try Gaussians from {C.SCALES[0]} to {C.SCALES[-1]:g} px and keep
the smallest one where the signal is significant above the background:</p>
<p class="eq">( G∗D − G∗B ) / √( max(G∗D, G∗B) · Σg² ) ≥ τ</p>
<p>Here D is the photon image, B is the background (rate × exposure), and Σg² turns the Poisson variance
into the variance of the smoothed value. Between two kernel sizes we interpolate, so no edges appear where
the size changes. Flux is G∗(D−B) / G∗E, with counts and exposure blurred the same way, which also takes
care of the mosaic edges. All blurs are done with FFTs, about 10 s per band for the full field.</p>
<p>On the first run the wind looked grainy. Some pixels passed the test at tiny kernels purely because of
upward noise fluctuations. A {C.SCALEMAP_MEDIAN}×{C.SCALEMAP_MEDIAN} median on the map of chosen sizes removes these single
pixels. Real sources cover several pixels, so they are not affected.</p>
<p>The brightness layer uses τ = {C.SNR_LUM:g} on the sum of all bands. The colour layer uses a stricter
τ = {C.SNR_COLOUR:g}, and the same kernel map is applied to every band, so the colour ratios always come from
equally smoothed data.</p>
<figure><img src="{F}fig3_smoothing.png"><figcaption>Fig. 3. The north-west part of the wind.
Fixed blurs are either noisy (b) or smear the sources (c). The adaptive version (d) is smooth in the gas
and sharp on the sources. (e) shows the kernel size used at each pixel.</figcaption></figure>
<table class="t"><tr><th>Method</th><th>Sky RMS</th><th>Wind S/N per pixel</th><th>Point source FWHM (px)</th></tr>{smooth_table}</table>
<p>The FWHM is measured on the source at {tuple(met['point_source_xy'])}, about 5′ off-axis, where the
telescope's point spread function is already a few pixels wide.</p>

<h2>6. Step 3: colour (<code>03_compose.py</code>)</h2>
<p>Brightness is the sum of all seven bands. Colour comes from three groups: 0.35–1.1 keV (the cooler wind
gas), 1.1–2.2 keV, and 2.2–6 keV (binaries and the absorbed nucleus). We show these as orange, green and
blue. We tried pure red, green and blue first, but the wind came out a flat, harsh red.</p>
<p>The colour layer is blurred another {C.COLOUR_BLUR:g} px. The eye notices colour detail much less than
brightness detail, the same idea as LRGB imaging. The whole galaxy is balanced to neutral, so a colour just
means softer or harder than the average. Following Lupton et al. (2004), only the colour ratios are applied
to the stretched brightness, so stretching does not shift the hues.</p>
<p><b>Reading the colours.</b> Orange is soft, thermal emission from oxygen and iron lines in the wind.
Toward the centre the colour moves through yellow and green to blue as the gas gets hotter and denser,
and as the dusty disk absorbs the softest photons. The blue haze along the disk is real hard emission from
unresolved X-ray binaries. The compact blue and white points are hard-spectrum binaries, including M82 X-1,
and background quasars seen through the galaxy.</p>
<figure><img src="{F}fig4_channels.png"><figcaption>Fig. 4. The three energy groups and the combined image.</figcaption></figure>

<h2>7. Step 4: stretch and finishing</h2>
<ol>
<li>Black point at {C.BLACK_NSIG}σ of the sky noise. White point at the {C.WHITE_PCT['wide']} percentile.</li>
<li>asinh stretch, softening {C.SOFT_NSIG['wide']}σ for the wide image and {C.SOFT_NSIG['core']}σ for the core, so both the faint gas and the
bright centre show.</li>
<li>Non-local-means denoising, applied only where the kernel was larger than {C.DETAIL_KERNEL:g} px, i.e. the gas and not the sources.</li>
<li>Local contrast (unsharp mask at σ = {C.LOCAL_CONTRAST[0]:g} px, amount {C.LOCAL_CONTRAST[1]}) to bring out the wind's filaments.</li>
<li>Sharpening (σ = {C.SHARPEN[0]} px, amount {C.SHARPEN[1]}), only where the kernel was below {C.DETAIL_KERNEL:g} px. Sharpening the noisy
areas would only have made up structure.</li>
<li>Saturation ×{C.SATURATION}, eased in with brightness so faint noisy areas stay neutral. The brightest parts fade
to white.</li>
<li>A 40 px fade at the mosaic edge, and export with north up as 16-bit TIFF, PNG and JPEG.</li>
</ol>
<figure><img src="{F}fig5_postprocess.png"><figcaption>Fig. 5. The core at each stage.</figcaption></figure>
<figure style="width: 60%; margin: 3mm auto;"><img src="{FI}M82_Chandra_core_Ramanujan.jpg">
<figcaption>Close-up of the centre (6.9′). The blue absorbed disk crosses the white core, with the wind
rising on both sides.</figcaption></figure>

<h2>8. Problems along the way</h2>
<table class="t"><tr><th>What we saw</th><th>Cause</th><th>Fix</th></tr>
<tr><td>Exposure map about 3× too high near the galaxy</td><td>Our first estimate used all bands, and the galaxy's soft halo leaked into the "background"</td><td>Used only the 2.8–6 keV band and masked a wider area</td></tr>
<tr><td>Streaks radiating from the centre of the exposure map</td><td>The way the masked hole was filled</td><td>Switched to a Laplace fill</td></tr>
<tr><td>A black spot, later a checkerboard, at X-1</td><td>The zero-count strip. Setting its exposure to zero made the flux unstable at small kernels</td><td>Filled the counts in from neighbouring pixels instead</td></tr>
<tr><td>Grainy wind</td><td>Noise pixels passing at tiny kernels, made worse by a minimum filter we had used</td><td>Median on the kernel map, τ raised from 4 to 6, denoising in the gas</td></tr>
<tr><td>Speckled colour and a purple haze</td><td>Noisy high-energy background clipped at zero</td><td>Stricter colour threshold, extra colour blur, less saturation in faint areas</td></tr>
</table>

<h2>9. Limitations</h2>
<ul>
<li>The exposure map is an estimate. The observation files needed to compute it properly were not provided.</li>
<li>We did not deconvolve. Sources far from the centre look larger because Chandra's resolution drops off-axis.</li>
<li>M82 X-1 is piled up, so its brightness and colour are not reliable.</li>
<li>X-rays have no colour. Ours are a mapping from energy, with orange for low and blue for high.</li>
</ul>

<h2>10. Code and reproduction</h2>
<p>All code is on GitHub at <a href="https://{REPO}">{REPO}</a>. The repository holds the five pipeline
scripts, the final images and figures, and a README with setup instructions. With the seven FITS files in
place, these commands rebuild every image, every figure and this PDF:</p>
<pre class="cmd">python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
./run_all.sh</pre>
<p>The run takes about two minutes. Every parameter is set in one file, <code>pipeline/config.py</code>, and all of
them are listed in the appendix. The software versions are in the next table.</p>
<table class="t"><tr><th>Software</th><th>Version</th><th>Use</th></tr>{versions}</table>

<h2>Credit</h2>
<p>X-ray data: NASA / CXC, Chandra X-ray Observatory (chandra.harvard.edu). The data were not captured by
Cepheid or by the authors. As the competition rules ask, this image must not be shared without this credit.</p>

<h2>References</h2>
<ul>
<li>Data: NASA / Chandra X-ray Center, chandra.harvard.edu. Not captured by Cepheid or by us.</li>
<li>Ebeling, White &amp; Rangarajan 2006, MNRAS 368, 65 (adaptive smoothing)</li>
<li>Lupton et al. 2004, PASP 116, 133 (colour composites)</li>
<li>Buades, Coll &amp; Morel 2005, CVPR (non-local means)</li>
<li>Strickland &amp; Heckman 2009, ApJ 697, 2030 (M82 superwind)</li>
<li>Lehnert, Heckman &amp; Weaver 1999, ApJ 523, 575 (the cap)</li>
<li>Bachetti et al. 2014, Nature 514, 202 (M82 X-2)</li>
</ul>

<h2 class="pb">Appendix: parameters</h2>
<table class="t"><tr><th>Name</th><th>Value</th></tr>{param_table}</table>
</body></html>"""

C.REPORT.mkdir(exist_ok=True)
HTML(string=html, base_url=str(C.REPORT)).write_pdf(C.REPORT / f"M82_Stacking_Report_{C.HOSTEL}.pdf")
