import numpy as np
from astropy.io import fits
from scipy import fft
from PIL import Image
import config as C


def band_name(b):
    return f"{b[0]}-{b[1]} eV"


def load_bands():
    files = [C.RAW_DIR / f"m82_{lo}_{hi}_eV.fits" for lo, hi in C.BANDS]
    cube = np.stack([fits.getdata(f).astype(np.float32) for f in files])
    return cube, fits.getheader(files[0])


def save_fits(name, data, hdr):
    fits.writeto(C.WORK / name, np.asarray(data, np.float32), hdr, overwrite=True)


def load_fits(name):
    return fits.getdata(C.WORK / name).astype(np.float32)


class FFTSmoother:
    """Gaussian blur at many sigmas from a single forward FFT (zero padded)."""

    def __init__(self, shape, pad=256):
        self.shape, self.pad = shape, pad
        self.py = fft.next_fast_len(shape[0] + 2 * pad)
        self.px = fft.next_fast_len(shape[1] + 2 * pad, real=True)
        fy = fft.fftfreq(self.py)[:, None]
        fx = fft.rfftfreq(self.px)[None, :]
        self.f2 = fy ** 2 + fx ** 2

    def forward(self, img):
        buf = np.zeros((self.py, self.px))
        p = self.pad
        buf[p:p + self.shape[0], p:p + self.shape[1]] = img
        return fft.rfft2(buf, workers=-1)

    def smooth(self, F, sigma):
        out = fft.irfft2(F * np.exp(-2 * np.pi ** 2 * sigma ** 2 * self.f2),
                         s=(self.py, self.px), workers=-1)
        p = self.pad
        return out[p:p + self.shape[0], p:p + self.shape[1]].astype(np.float32)


def kernel_sum_sq(sigma):
    # sum of g^2 for a normalised discrete gaussian -> Poisson variance factor
    x = np.arange(-int(6 * sigma) - 2, int(6 * sigma) + 3)
    g = np.exp(-0.5 * (x / sigma) ** 2)
    g /= g.sum()
    return float((g ** 2).sum() ** 2)


def asinh_stretch(x, soft, vmax):
    return np.arcsinh(np.clip(x, 0, None) / soft) / np.arcsinh(vmax / soft)


def to_uint(img, bits=8):
    top = 255 if bits == 8 else 65535
    return (np.clip(img, 0, 1) * top + 0.5).astype(np.uint8 if bits == 8 else np.uint16)


def save_png(path, img):
    # row 0 of the array is south, so flip for north up
    Image.fromarray(to_uint(img[::-1])).save(path)
