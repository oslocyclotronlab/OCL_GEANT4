import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
from scipy.special import erfc
from tabulate import tabulate
from typing import Dict, Tuple, Optional, Callable

import fnmatch
import os

# from uncertainties import ufloat
from uncertainties import unumpy, ufloat

from ompy import Matrix


class Recalibrate:
    def __init__(self):
        self._popt = self.get_coefficients()
        pass

    def __call__(self, x):
        return self.func(x, *self._popt)

    @staticmethod
    def func(x, *p):
        pol = np.poly1d(p[0:3])
        x_less = pol(x)
        pol = np.poly1d(p[3:5])
        x_larger = pol(x)
        return np.where(x<500, x_less, x_larger)

    @staticmethod
    def get_coefficients(plot=False):
        data = np.loadtxt("coords_Eu.txt")
        x1 = data[1::2, 0]
        y1 = data[0::2, 0]
        if plot:
            fig, ax = plt.subplots()
            ax.plot(x1, y1, "o")

        data = np.loadtxt("coords_Co.txt")
        x2 = data[1::2, 0]
        y2 = data[0::2, 0]
        if plot:
            ax.plot(x2, y2, "o")

        x = np.append(x1, x2)
        y = np.append(y1, y2)

        x = np.sort(x)
        y = np.sort(y)

        res = np.polyfit(x, y, 1)
        p = np.poly1d(res)
        # print(p)
        if plot:
            ax.plot(x, p(x), "-", alpha=0.2)

        res = np.polyfit(x, y, 2)
        p = np.poly1d(res)
        # print(p)
        if plot:
            ax.plot(x, p(x), "--", alpha=0.2)

        popt, pcov = curve_fit(Recalibrate.func, x, y, p0=[-7e7, 1, -1, 1, 1])
        # def func(x, *p):
        #     pol = np.poly1d(p[0:3])
        #     x_less = pol(x) + p[3]*np.sqrt(x)
        #     pol = np.poly1d(p[4:6])
        #     x_larger = pol(x)
        #     return np.where(x<500, x_less, x_larger)

        # popt, pcov = curve_fit(func, x, y, p0=[ -7e7, 1, -1, 1e-1, 1, 1])

        # print(popt)
        if plot:
            ax.plot(x, Recalibrate.func(x, *popt), "-", alpha=0.2)

        return popt


calibrate = Recalibrate()


class SpectrumComparison:

    def __init__(self):
        pass

    def get_data(self, fname_sim, fname_exp, fname_bg, fwhm_pars,
                 measure_time_exp, measure_time_bg, idet=0,
                 recalibrate=True):

        sim_all = Matrix(path=fname_sim)
        sim = np.zeros((len(sim_all.Eg), 2))
        sim[:, 0] = sim_all.Eg
        sim[:, 0] *= 1000  # MeV to keV
        if isinstance(idet, int):
            sim[:, 1] = sim_all.values[idet, :]
        else:
            sim[:, 1] = sim_all.values[idet, :].sum(axis=0)
        sim = do_smooth(sim, fwhm_pars)
        self.fwhm_pars = fwhm_pars

        # exp_and_bg = np.loadtxt(fname_exp)
        bg_all = Matrix(path=fname_bg)
        bg = np.zeros((len(bg_all.Eg), 2))
        bg[:, 0] = bg_all.Eg
        if isinstance(idet, int):
            bg[:, 1] = bg_all.values[idet, :]
        else:
            bg[:, 1] = bg_all.values[idet, :].sum(axis=0)

        # subtract bgfrom experiment
        exp_all = Matrix(path=fname_exp)
        exp = np.zeros((len(exp_all.Eg), 2))
        exp[:, 0] = exp_all.Eg
        xexp = exp_all.Eg
        if isinstance(idet, int):
            yexp = exp_all.values[idet, :]
        else:
            yexp = exp_all.values[idet, :].sum(axis=0)

        if recalibrate:
            xexp = calibrate(xexp)
            exp[:, 0] = xexp

        # exp = np.copy(exp_and_bg)
        # xexp = exp[:, 0]
        # yexp = exp[:, 1]
        uyexp = unumpy.uarray(yexp, np.sqrt(yexp))
        uybg = unumpy.uarray(bg[:, 1], np.sqrt(bg[:, 1]))

        bg_ratio = measure_time_exp / measure_time_bg
        uyexp = uyexp - uybg * bg_ratio
        exp[:, 1] = unumpy.nominal_values(uyexp)

        # fig, ax = plt.subplots()
        # ax.semilogy(xexp, exp[:, 1]/70)
        # ax.plot(sim[:, 0], sim[:, 1])

        self.sim = sim
        self.exp = exp
        self.xexp = xexp
        self.xsim = sim[:, 0]
        self.uyexp = uyexp

    def scale_sim_to_exp_fit(self, E1, E2):
        exp = self.exp
        sim = self.sim
        (popt_exp, pcov_exp), xfit = do_fit(E1, E2, exp,
                                            fwhm_pars=self.fwhm_pars)
        (popt_sim, pcov_sim), xfit = do_fit(E1, E2, sim,
                                            fwhm_pars=self.fwhm_pars)

        ratio = popt_exp / popt_sim
        ratio = ratio[0]
        # print(popt_exp[0])
        # alternative: scale the number of counts (above some threshold)
        # print("ratio; scaling all counts", ratio, exp[100:, 1].sum()/sim[20:, 1].sum()/5)
        # ratio = exp[100:, 1].sum()/sim[20:, 1].sum()/5

        # print(sim)
        # xsim = sim[:, 0]
        ysim = sim[:, 1]
        uysim_scaled = unumpy.uarray(ysim, np.sqrt(ysim))
        uysim_scaled *= ratio

        self.uysim_scaled = uysim_scaled
        self.fsim_scaled = uinterp1D(self.xsim, uysim_scaled)
        self.fexp = uinterp1D(self.xexp, self.uyexp)

        self.xfit = {"peak_fit": xfit}
        self.scale_factor = ratio
        self.popt_sim = popt_sim
        self.popt_exp = popt_exp

    def scale_sim_to_exp_area(self, E1, E2):
        exp = self.exp
        sim = self.sim

        counts_exp = self.get_area(exp, E1, E2)
        counts_sim = self.get_area(sim, E1, E2)

        iE1 = np.abs(sim[:, 0]-E1).argmin()
        iE2 = np.abs(sim[:, 0]-E2).argmin()
        counts_sim = sim[iE1:iE2+1, 1].sum()

        binwidth_sim = sim[1, 0] - sim[0, 0]
        binwidth_exp = exp[1, 0] - exp[0, 0]
        ratio = counts_exp / counts_sim * (binwidth_exp/binwidth_sim)

        # print(sim)
        # xsim = sim[:, 0]
        ysim = sim[:, 1]
        uysim_scaled = unumpy.uarray(ysim, np.sqrt(ysim))
        uysim_scaled *= ratio

        self.uysim_scaled = uysim_scaled
        self.fsim_scaled = uinterp1D(self.xsim, uysim_scaled)
        self.fexp = uinterp1D(self.xexp, self.uyexp)

        self.scale_factor = ratio
        self.xfit = {"area": [E1, E2]}

    @staticmethod
    def get_area(arr, E1, E2):
        iE1 = np.abs(arr[:, 0]-E1).argmin()
        iE2 = np.abs(arr[:, 0]-E2).argmin()
        return arr[iE1:iE2+1, 1].sum()


    def scale_sim_to_exp_manual(self, ratio):
        exp = self.exp
        sim = self.sim

        # print(sim)
        # xsim = sim[:, 0]
        ysim = sim[:, 1]
        uysim_scaled = unumpy.uarray(ysim, np.sqrt(ysim))
        uysim_scaled *= ratio

        self.uysim_scaled = uysim_scaled
        self.fsim_scaled = uinterp1D(self.xsim, uysim_scaled)
        self.fexp = uinterp1D(self.xexp, self.uyexp)

        self.scale_factor = ratio
        self.xfit = {"ratio": ratio}

    def get_chi2(self):
        xexp = self.xexp
        fexp = self.fexp
        fsim_scaled = self.fsim_scaled

        yexp = fexp(xexp)
        # yexp[yexp == 0] = np.nan

        chi2 = (unumpy.nominal_values(yexp) -
                unumpy.nominal_values(fsim_scaled(xexp)))**2
        chi2 /= unumpy.std_devs(yexp)**2 \
            + unumpy.std_devs(fsim_scaled(xexp))**2
        self.chi2 = chi2
        return chi2

    def get_rel_diff(self, smooth_window_keV=20):
        xexp = self.xexp
        fexp = self.fexp
        fsim_scaled = self.fsim_scaled

        yexp = fexp(xexp)
        yexp[unumpy.nominal_values(yexp) == 0] = ufloat(np.nan, np.nan)

        rel_diff = (yexp - fsim_scaled(xexp)) / yexp * 100
        rel_diff_err = unumpy.std_devs(rel_diff)
        rel_diff = unumpy.nominal_values(rel_diff)

        rel_diff_smooth = unumpy.nominal_values(yexp - fsim_scaled(xexp))
        window_len = smooth_window_keV / (xexp[1] - xexp[0])
        rel_diff_smooth = smooth(rel_diff_smooth, window_len=int(window_len-1))
        rel_diff_smooth /= unumpy.nominal_values(yexp) / 100

        self.rel_diff = rel_diff
        self.rel_diff_err = rel_diff_err
        self.rel_diff_smooth = rel_diff_smooth
        return rel_diff, rel_diff_smooth

    def plots(self, xmax=None, ymin=1.2e2, ymax=1e6, title=None,
              plot_smoothed=True):

        xexp = self.xexp

        if hasattr(self, 'rel_diff'):
            pass
        else:
            self.get_rel_diff()

        # fig, (ax, ax2) = plt.subplots(2, 1, sharex=True,
        #                               gridspec_kw={'hspace': 0.0})
        fig = plt.figure(figsize=[6.4, 3.5], dpi=200)
        ax = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
        ax2 = plt.subplot2grid((3, 1), (2, 0), sharex=ax)
        fig.subplots_adjust(hspace=0, bottom=0.145)
        plt.setp(ax.get_xticklabels(), visible=False)

        ax.semilogy(xexp, self.exp[:, 1], label="exp.")
        ax.plot(self.xsim, unumpy.nominal_values(self.uysim_scaled),
                label="sim., scaled")

        if "peak_fit" in self.xfit:
            xfit = self.xfit["peak_fit"]
            ax.plot(xfit, gaus_with_expbg(xfit, *self.popt_exp),
                    label="fit_to_exp")
            ax.plot(xfit,
                    gaus_with_expbg(xfit, *self.popt_sim) * self.scale_factor,
                    label="fit_to_sim")
        elif "area" in self.xfit:
            E1, E2 = self.xfit["area"]
            ax.axvspan(E1, E2, color="k", alpha=0.2, lw=0,
                       label="scaling area")
        elif "ratio" in self.xfit:
            pass

        # y0 = unumpy.nominal_values(self.rel_diff)
        # yerr = unumpy.std_devs(self.rel_diff)
        y0 = self.rel_diff
        yerr = self.rel_diff_err
        ax2.fill_between(xexp, y0 - yerr, y0 + yerr, color="k", alpha=0.2)
        ax2.plot(xexp, y0, c="k", alpha=0.6)
        if plot_smoothed:
            ax2.plot(xexp, self.rel_diff_smooth, "C5--", alpha=0.8,
                     label="smoothed")

        ax2.axhline(y=0, color="r", lw=1)
        # ax2.set_yscale('symlog')
        # ax2.set_xlim(0, 1400)

        ax.set_yscale("log")
        ax.set_xlabel("Energy [keV]")
        ax.set_ylabel("counts / bin")
        # ax.legend(loc="best")
        ax.set_xlim(0, xmax)
        ax.set_ylim(ymin, ymax)

        # ax2.legend(loc="best")
        ax2.set_ylim(-18, 18)
        ax2.set_xlabel("Energy [keV]")
        ax2.set_ylabel(r"$\frac{\mathrm{exp.}-\mathrm{sim.}}"
                       "{\mathrm{exp.}} [\%]$")

        ax2.tick_params(labelbottom='on', top='on')

        if title is not None:
            fig.suptitle(title, fontsize=12)
        return fig, (ax, ax2)

    @staticmethod
    def cut_array(xarr, yarr, Elow, Ehigh):
        assert(Elow < Ehigh), "Elow has to be larger then Ehigh"
        idx1 = np.abs(xarr - Elow).argmin()
        idx2 = np.abs(xarr - Ehigh).argmin()
        return yarr[idx1:idx2 + 1]

    def fom(self, Ecompare_low, Ecompare_high, printout=True):
        if hasattr(self, 'rel_diff'):
            pass
        else:
            self.get_rel_diff()

        self.get_chi2()
        chi2 = np.copy(self.chi2)
        rel_diff = np.copy(self.chi2)
        rel_diff_smooth = np.copy(self.chi2)
        chi2 = self.cut_array(self.xexp, chi2, Ecompare_low, Ecompare_high)
        chi2 = chi2.sum()
        rel_diff = self.cut_array(self.xexp, self.rel_diff, Ecompare_low,
                                  Ecompare_high)
        rel_diff = abs(rel_diff).mean()
        rel_diff_smooth = self.cut_array(self.xexp, self.rel_diff_smooth,
                                         Ecompare_low, Ecompare_high)

        rel_diff_smooth = abs(rel_diff_smooth).mean()

        if printout:
            print("chi2 between {} and {} keV: {}".
                  format(Ecompare_low, Ecompare_high, chi2))
            print("Mean abs(rel_diff) between {} and {} keV: {:.2f}%"
                  .format(Ecompare_low, Ecompare_high, rel_diff))
            print("Mean smoothed abs(rel_diff) between {} and {} keV: {:.2f}%"
                  .format(Ecompare_low, Ecompare_high, rel_diff_smooth))

        return chi2, rel_diff, rel_diff_smooth


def fFWHM(E, p):
    # sgima in keV
    return np.sqrt(p[0] + p[1] * E + p[2] * E**2)


def fsigma(E, p):
    # sgima in keV
    return fFWHM(E, p) / (2 * np.sqrt(2 * np.log(2)))


def do_smooth(data, pars):
    data = np.copy(data)
    Ebins = data[:, 0]
    counts = data[:, 1]
    nbins = len(data)
    smoothed = np.zeros(nbins)
    binwidth = Ebins[1] - Ebins[0]
    for i, E in enumerate(Ebins):
        if E > 0:
            # FWHM=fFHWM(E, pars)/binwidth
            # sigma = np.sqrt(pars[0] + pars[1]*E + pars[2]*E**2)
            # sigma = FWHM/(2*np.sqrt(2*np.log(2)))
            sigma = fsigma(E, pars) / binwidth
            # smooth each bin in the spectrum: create new "spectrum",
            # which is filled only at E
            if counts[i] != 0:  # smooth this!
                counts_ = np.zeros(nbins)
                counts_[i] = counts[i]
                smoothed += gaussian_filter1d(counts_, sigma)
    data[:, 1] = smoothed
    return data


def exp_scaled(x, scale, decay_const):
    return scale * np.exp(-decay_const * x)


def gauss(x, scale, loc, sigma):
    return scale / (2. * np.pi * sigma) * np.exp(-(x - loc)**2 / (2 * sigma**2))


def gaus_with_expbg(x, scale, loc, sigma, bg_scale, bg_decay):
    line = gauss(x, scale, loc, sigma)
    bg = exp_scaled(x, bg_scale, bg_decay)
    return line + bg


def guess_p0(data, fitfun=gaus_with_expbg, fwhm_pars=None):
    assert(fitfun == gaus_with_expbg), "Other Bg estimation not implemented"

    data = np.copy(data)
    # gaus_with_expbg
    # scale, loc, sigma, bg_scale, bg_decay
    x1, y1 = data[0, :]
    x2, y2 = data[-1, :]
    if y1 < 1e-4:  # workaround for exp data with too much bg
        y1 = 1e-4
    if y2 < 1e-4:  # workaround for exp data with too much bg
        y2 = 1e-4
    bg_decay = (np.log(y2) - np.log(y1)) / (x2 - x1)
    bg_decay = -bg_decay
    bg_scale = y1 / np.exp(-bg_decay * x1)

    bg = exp_scaled(data[:, 0], bg_scale, bg_decay)
    data[:, 1] -= bg

    scale = data[:, 1].max()
    loc = data[data[:, 1] == scale].flatten()[0]
    sigma = fsigma(loc, fwhm_pars)

    # peak cross section to area
    scale *= (2. * np.pi * sigma)
    return [scale, loc, sigma, bg_scale, bg_decay]


def do_fit(Elow, Ehigh, data,
           p0=None,
           fitfunc=gaus_with_expbg, bounds="positive",
           fwhm_pars=None):
    idx1 = np.abs(data[:, 0] - Elow).argmin()
    idx2 = np.abs(data[:, 0] - Ehigh).argmin()
    xdata = data[idx1:idx2 + 1, 0]
    ydata = data[idx1:idx2 + 1, 1]

    if p0 is None:
        p0 = guess_p0(np.c_[xdata, ydata], fwhm_pars=fwhm_pars)
        # plt.plot(xdata, ydata)
        # plt.plot(xdata, fitfunc(xdata, *p0))
        # plt.show()

    if bounds == "positive":
        p_up = np.copy(p0)
        p_up[:] = np.inf
        bounds = (np.zeros(len(p0)), p_up)
        bounds[0][1] = Elow
        bounds[1][1] = Ehigh
    else:
        bounds = bounds

    # convert sigma to bin-width
    # binwidth = data[1, 0] - data[0, 0]
    # p0[2] /= binwidth
    return curve_fit(fitfunc, xdata, ydata, p0=p0, bounds=bounds), xdata


def uinterp1D(x, uvar, fill_value="extrapolate"):
    fnom = interp1d(x, unumpy.nominal_values(uvar), fill_value=fill_value)
    fstd = interp1d(x, unumpy.std_devs(uvar), fill_value=fill_value)

    print()
    def f(var):
        std = fstd(var)
        std[std < 0] = np.nan
        # std = fstd(var)[fstd(var) < 0] = np.nan
        return unumpy.uarray(fnom(var), std)
    return f


def smooth(x, window_len, window='hanning', trim=True):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window;
        should be an odd integer
        : the type of window from 'flat', 'hanning',
        'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
        trim: trim output so same length as imput

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett,
    numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead
     of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError(
            "Window is on of 'flat', 'hanning', 'hamming', 'bartlett',"
            " 'blackman'")

    s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')
    if trim:
        if int(window_len) % 2:  # uneven
            upper = -int(window_len / 2)-1
        else:
            upper = -int(window_len / 2)
        return y[int(window_len / 2 - 1):upper]
    else:
        return y


class PeakFitter:
    def __init__(self, x: np.ndarray, y: np.ndarray,
                 yerr: Optional[np.ndarray] = None):
        self.x = x.copy()
        self.y = y.copy()

        if yerr is None:
            self.yerr = None
        else:
            self.yerr = yerr.copy()
            self.yerr[self.y == 0] = np.inf

        self._xfit = None  # make it easier to retrieve current fit x values

    def fitGausBg(self, pre_region: Tuple[float, float],
                  post_region: Tuple[float, float]) -> Dict[str, float]:
        """Fit a gamma peak with a gaussian on top of a const. bg

        Args:
            pre_region: Energy region with linear spectra on the left
            side of the peak
            post_region: Energy region with linear spectra on the right
            side of the peak

        Returns:
            popt, cov: optimal parameters and covariance
        """
        x = self.x
        y = self.y

        pre_slice, peak_slice, post_slice, fit_slice = \
            self.get_fit_slices(pre_region, post_region)

        # pre_spec = y[pre_slice]
        peak_spec = y[peak_slice]
        # post_spec = y[post_slice]

        # Estimate the mean, std and constant
        peak_mean = np.sum(x[peak_slice]*peak_spec)/np.sum(peak_spec)
        peak_var = np.sum(peak_spec*x[peak_slice]**2)/np.sum(peak_spec)
        peak_std = np.sqrt(peak_var - peak_mean**2)

        # We estimate the gauss constant from the height found at the mean
        peak_const = np.max(peak_spec) / self.gaus(peak_mean, 1., peak_mean,
                                                   peak_std)

        # pol_estimate = np.polyfit(np.append(x[pre_slice], x[post_slice]),
        #                           np.append(y[pre_slice], y[post_slice]), 0)

        # Calculate the curve fitting
        initial_guess = [peak_const, peak_mean, peak_std,
                         y[post_slice].mean()]
        yerr = self.yerr[fit_slice] if self.yerr is not None else None
        popt, cov = curve_fit(self.gaus_bg, x[fit_slice], y[fit_slice],
                              p0=initial_guess, sigma=yerr)
        return popt, cov

    def fitDoubleGausBg(self, pre_region: Tuple[float, float],
                        post_region: Tuple[float, float],
                        p0_E: Tuple[float, float],
                        p0_sigma=np.array([10., 10.]),
                        ) -> Dict[str, float]:
        """Fit a two gaussians on top of a const. bg

        Args:
            pre_region: Energy region with linear spectra on the left
            side of the peak
            post_region: Energy region with linear spectra on the right
            side of the peak
            p0_E: Initial guess for peaks
            p0_sigma: Initial guess for sigma

        Returns:
            popt, cov: optimal parameters and covariance
        """
        x = self.x
        y = self.y

        pre_slice, peak_slice, post_slice, fit_slice = \
            self.get_fit_slices(pre_region, post_region)

        peak_spec = y[peak_slice]

        # We estimate the constant from the height found at the mean
        p0_const = np.max(peak_spec) / self.gaus(p0_E[0], 1., p0_E[0],
                                                 p0_sigma[0])
        p0_const2 = np.max(peak_spec) / self.gaus(p0_E[1], 1., p0_E[1],
                                                  p0_sigma[1])

        # pol_estimate = np.polyfit(np.append(x[pre_slice], x[post_slice]),
        #                           np.append(y[pre_slice], y[post_slice]), 0)

        # Calculate the curve fitting
        initial_guess = [p0_const, p0_E[0], p0_sigma[0],
                         p0_const2, p0_E[1], p0_sigma[1],
                         y[post_slice].mean()]
        yerr = self.yerr[fit_slice] if self.yerr is not None else None
        popt, cov = curve_fit(self.doublegaus_bg, x[fit_slice], y[fit_slice],
                              p0=initial_guess, sigma=yerr)
        return popt, cov

    def fitGausBgStep(self, pre_region: Tuple[float, float],
                      post_region: Tuple[float, float]) -> Dict[str, float]:
        """Fit a gamma peak with a gaussian on top of a const. bg + step fct

        Args:
            pre_region: Energy region with linear spectra on the left
            side of the peak
            post_region: Energy region with linear spectra on the right
            side of the peak
        """
        x = self.x
        y = self.y

        pre_slice, peak_slice, post_slice, fit_slice = \
            self.get_fit_slices(pre_region, post_region)

        popt_gaus, _ = self.fitGausBg(pre_region, post_region)
        # Calculate the curve fitting
        p0_step_constant = 1000
        initial_guess = [*popt_gaus, p0_step_constant]
        yerr = self.yerr[fit_slice] if self.yerr is not None else None
        popt, cov = curve_fit(self.gaus_bg_step, x[fit_slice], y[fit_slice],
                              p0=initial_guess, sigma=yerr,
                              # bounds=[0, np.inf]
                              )
        return popt, cov

    def fitGausStep(self, pre_region: Tuple[float, float],
                    post_region: Tuple[float, float]) -> Dict[str, float]:
        """Fit a gamma peak with a gaussian on top and step fct

        Args:
            pre_region: Energy region with linear spectra on the left
            side of the peak
            post_region: Energy region with linear spectra on the right
            side of the peak
        """
        x = self.x
        y = self.y

        pre_slice, peak_slice, post_slice, fit_slice = \
            self.get_fit_slices(pre_region, post_region)

        popt_gaus, _ = self.fitGausBg(pre_region, post_region)
        # Calculate the curve fitting
        p0_step_constant = 1000
        initial_guess = [*popt_gaus[:-1], p0_step_constant]
        yerr = self.yerr[fit_slice] if self.yerr is not None else None
        popt, cov = curve_fit(self.gaus_step, x[fit_slice], y[fit_slice],
                              p0=initial_guess, sigma=yerr,
                              # bounds=[0, np.inf]
                              )
        return popt, cov

    def fitDoubleGausBgStep(self, pre_region: Tuple[float, float],
                            post_region: Tuple[float, float],
                            p0_E: Tuple[float, float],
                            p0_sigma=np.array([10., 10.]),
                            ) -> Dict[str, float]:
        """Fit two gaussians on top of a const. bg + step fct

        Args:
            pre_region: Energy region with linear spectra on the left
            side of the peak
            post_region: Energy region with linear spectra on the right
            side of the peak
            p0_E: Initial guess for peaks
            p0_sigma: Initial guess for sigma
        """
        x = self.x
        y = self.y

        pre_slice, peak_slice, post_slice, fit_slice = \
            self.get_fit_slices(pre_region, post_region)

        popt_gaus, _ = self.fitDoubleGausBg(pre_region, post_region,
                                            p0_E, p0_sigma)
        # Calculate the curve fitting
        p0_step_constant = 0.1
        initial_guess = [*popt_gaus, p0_step_constant]
        yerr = self.yerr[fit_slice] if self.yerr is not None else None

        bounds = [[0, np.inf] for i in initial_guess]

        popt, cov = curve_fit(self.doublegaus_bg_step,
                              x[fit_slice], y[fit_slice],
                              p0=initial_guess, sigma=yerr,
                              # bounds=[0, np.inf]
                              )
        return popt, cov

    def fitGausLinBgStep(self, pre_region: Tuple[float, float],
                         post_region: Tuple[float, float]) -> Dict[str, float]:
        """Fit a gamma peak with a gaussian on top of a linear bg + step fct

        Args:
            pre_region: Energy region with linear spectra on the left
            side of the peak
            post_region: Energy region with linear spectra on the right
            side of the peak
        """
        x = self.x
        y = self.y

        pre_slice, peak_slice, post_slice, fit_slice = \
            self.get_fit_slices(pre_region, post_region)

        popt_gaus, _ = self.fitGausBgStep(pre_region, post_region)
        # Calculate the curve fitting
        pol_estimate = np.polyfit(np.append(x[pre_slice], x[post_slice]),
                                  np.append(y[pre_slice], y[post_slice]), 1)[0]

        initial_guess = [*popt_gaus[:-1], pol_estimate, popt_gaus[-1]]
        yerr = self.yerr[fit_slice] if self.yerr is not None else None
        popt, cov = curve_fit(self.gaus_linbg_step, x[fit_slice], y[fit_slice],
                              p0=initial_guess, sigma=yerr,
                              # bounds=[0, np.inf]
                              )
        return popt, cov

    def get_fit_slices(self, pre_region: Tuple[float, float],
                       post_region: Tuple[float, float]):
        """ Extract fit slices from tuple of Emin/Emac before / after peak

        Args:
            pre_region: Energy region with linear spectra on the left
            side of the peak
            post_region: Energy region with linear spectra on the right
            side of the peak

        Returns:
            pre_slice, peak_slice, post_slice, fit_slice
        """
        x = self.x

        def ix(x0):
            return (np.abs(x-x0)).argmin()

        pre_slice = slice(ix(pre_region[0]), ix(pre_region[1])+1)
        peak_slice = slice(ix(pre_region[1])+1, ix(post_region[0]))
        post_slice = slice(ix(post_region[0]), ix(post_region[1])+1)
        fit_slice = slice(ix(pre_region[0]), ix(post_region[1])+1)

        self._xfit = x[fit_slice]
        return pre_slice, peak_slice, post_slice, fit_slice

    @staticmethod
    def gaus(x, const, mean, std):
        """ Gaussian """
        return const/(std*np.sqrt(2*np.pi)) * np.exp(-0.5*((x-mean)/std)**2)

    @staticmethod
    def smoothstep(x, const, mean, std):
        """ Smooth step """
        return const * erfc((x - mean) / (np.sqrt(2)*std))

    @staticmethod
    def gaus_bg(x, gaus_const, mean, std, offsett):
        """ Gauss and constant background """
        return PeakFitter.gaus(x, gaus_const, mean, std) + offsett

    @staticmethod
    def doublegaus_bg(x, gaus_const, mean, std,
                      gaus_const2, mean2, std2, offsett):
        """ Gauss and constant background """
        return (PeakFitter.gaus(x, gaus_const, mean, std)
                + PeakFitter.gaus(x, gaus_const2, mean2, std2) + offsett)

    @staticmethod
    def gaus_linbg(x, gaus_const, mean, std, offsett, bgslope):
        """ Gauss and linear background """
        return PeakFitter.gaus(x, gaus_const, mean, std) + x*bgslope + offsett

    @staticmethod
    def gaus_bg_step(x, gaus_const, mean, std, offsett, step_const):
        """ Gauss and constant background and step function """
        gaus_bg = PeakFitter.gaus_bg(x, gaus_const, mean, std, offsett)
        step = PeakFitter.smoothstep(x, step_const, mean, std)
        return gaus_bg + step

    @staticmethod
    def gaus_step(x, gaus_const, mean, std, step_const):
        """ Gauss and step function """
        gaus = PeakFitter.gaus(x, gaus_const, mean, std)
        step = PeakFitter.smoothstep(x, step_const, mean, std)
        return gaus + step

    @staticmethod
    def doublegaus_bg_step(x, gaus_const, mean, std, gaus_const2, mean2, std2,
                           offsett, step_const):
        """ Gauss and constant background and step function """
        doublegaus_bg = PeakFitter.doublegaus_bg(x, gaus_const, mean, std,
                                                 gaus_const2, mean2, std2,
                                                 offsett)
        step = PeakFitter.smoothstep(x, step_const, mean, std)
        return doublegaus_bg + step

    @staticmethod
    def gaus_linbg_step(x, gaus_const, mean, std, offsett, bgslope,
                        step_const):
        """ Gauss and linaer background and step function """
        gaus_linbg = PeakFitter.gaus_linbg(x, gaus_const, mean, std, offsett,
                                           bgslope)
        step = PeakFitter.smoothstep(x, step_const, mean, std)
        return gaus_linbg + step
