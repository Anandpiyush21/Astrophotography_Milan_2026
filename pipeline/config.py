from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "Astrophotography Raw FITS"
WORK = ROOT / "work"
OUT = ROOT / "outputs"
REPORT = ROOT / "report"

HOSTEL = "Ramanujan"
TEAM = ["Piyush Anand", "Darshanraj Pattanaik"]

BANDS = [(350, 500), (500, 700), (700, 1100), (1100, 1600),
         (1600, 2200), (2200, 2800), (2800, 6000)]  # eV
PIXEL_ARCSEC = 0.492
DISTANCE_MPC = 3.6

# exposure / background
EXPOSURE_BAND = 6        # 2.8-6 keV, mostly particle background
BKG_BIN = 16
BKG_MEDIAN_BOX = 5
GALAXY_MASK_R = 480      # px
BKG_FIT_R = 1200         # px, fit background only outside this radius
MIN_REL_EXPOSURE = 0.2

# adaptive smoothing
SCALES = [0.5, 0.7, 1.0, 1.4, 2.0, 2.8, 4.0, 5.6, 8.0, 11.0, 16.0, 22.0, 32.0]
SNR_LUM = 6.0
SNR_COLOUR = 15.0
SCALEMAP_MEDIAN = 5
SCALEMAP_SMOOTH = 1.5

# colour
RGB_GROUPS = {"R": [0, 1, 2], "G": [3, 4], "B": [5, 6]}
GROUP_COLOURS = {"R": (1.00, 0.42, 0.10), "G": (0.30, 1.00, 0.35), "B": (0.15, 0.40, 1.00)}
COLOUR_BLUR = 5.0
SATURATION = 1.45
SAT_RAMP = (0.08, 0.35)
HIGHLIGHT_START = 0.82

# stretch and finishing
WIDE_BOX = (540, 2740, 620, 2820)   # x0, x1, y0, y1
CORE_HALF = 420
BLACK_NSIG = 0.5
SOFT_NSIG = {"wide": 1.1, "core": 12.0}
WHITE_PCT = {"wide": 99.995, "core": 99.98}
NLM = {"patch_size": 5, "patch_distance": 7, "h_factor": 0.9}
LOCAL_CONTRAST = (40.0, 0.35)       # sigma, amount
SHARPEN = (1.2, 0.40)               # sigma, amount
DETAIL_KERNEL = 2.0                 # px, only sharpen where kernel is smaller than this
