import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
import fnmatch
import os
import re
from tqdm import tqdm
from datetime import datetime

import sys
from uncertainties import unumpy

import pandas as pd
from ompy import Matrix, Vector
# sys.path.append("../exp")

from compare import SpectrumComparison, do_smooth


class SpectrumComparisonInBeam(SpectrumComparison):
    def __init__(self):
        pass

    def get_data(self, exp, bg, bg_ratio, sim, fwhm_pars):
        self.fwhm_pars = fwhm_pars

        xexp = exp[:, 0]
        yexp = exp[:, 1]
        uyexp = unumpy.uarray(yexp, np.sqrt(yexp))
        uybg = unumpy.uarray(bg[:, 1], np.sqrt(bg[:, 1]))

        uyexp = uyexp - uybg * bg_ratio
        exp[:, 1] = unumpy.nominal_values(uyexp)

        sim = do_smooth(sim, fwhm_pars)

        self.sim = sim
        self.exp = exp
        self.xexp = exp[:, 0]
        self.xsim = sim[:, 0]
        self.uyexp = uyexp

    def plot_current(self):
        fig, ax = plt.subplots()
        ax.semilogy(self.exp[:, 0], self.exp[:, 1])
        ax.plot(self.sim[:, 0], self.sim[:, 1])
        plt.show()


def get_fom(fnamesim,
            exp, bg, fwhm_pars,
            bg_ratio,
            Efit_low, Efit_high,
            do_plot=True, printout=False,
            manual_ratio=None, xmax_plot=1400,
            ymin=1.2e3, ymax=1e7):
    """ get figure of merrit

    fnamesim: str, like "60Co", or "152Eu", "1777keV"
    """
    fname_sims = []
    for file in os.listdir('mama_spectra/root_files'):
        if fnmatch.fnmatch(file, f'inbeam_grid_*{fnamesim}*.m'):
            if "_all.m" in file:
                continue
            fname_sims.append(os.path.join("mama_spectra/root_files", file))
    fname_sims.sort()
    grid_points = np.full_like(fname_sims, np.nan, dtype=float)
    foms = np.zeros((len(fname_sims), 3))

    for i, fname_sim in enumerate(tqdm(fname_sims)):
        # if i > 2:
        #     break
        if printout:
            print("fitting: ", fname_sim)
        sc = SpectrumComparisonInBeam()

        sim = Vector(path=fname_sim)
        sim = np.c_[sim.E*1000, sim.values]

        sc.get_data(exp.copy(), bg.copy(), bg_ratio, sim, fwhm_pars)

        if manual_ratio is None:
            sc.scale_sim_to_exp_area(Efit_low, Efit_high)
        else:
            sc.scale_sim_to_exp_manual(manual_ratio)

        # sc.plot_current()
        chi2 = sc.get_chi2()
        rel_diff, rel_diff_smooth = sc.get_rel_diff(smooth_window_keV=20)

        foms[i, :] = sc.fom(Ecompare_low, Ecompare_high, printout=False)
        # print(sc.fom(Ecompare_low, Ecompare_high))
        grid_points[i] = int(re.search(r"grid_(-*\d*)_", fname_sim)[1])

        if do_plot:
            fig, (ax1, ax2) = sc.plots(title=fname_sim, xmax=xmax_plot,
                                       plot_smoothed=False)
        ax1.set_ylim(ymin, ymax)
        fig.savefig(f"figs_inbeam/{fnamesim}_{grid_points[i]:.0f}_noleg.png")
        ax1.legend()
        fig.savefig(f"figs_inbeam/{fnamesim}_{grid_points[i]:.0f}.png")
        # plt.show()
        plt.close(fig)


def arr_from_py(fname, Emin, Emax, remove_negative=False):
    mat = Matrix(path=fname)
    # mat.plot(vmin=1, vmax=1e4)
    if remove_negative:
        mat.remove_negative()
    values, E = mat.projection("Eg", Emin, Emax)

    # mat1 = mat.copy()
    # mat2 = mat.copy()
    # mat1.cut("Eg", 550, 550)
    # mat2.cut("Eg", 1778-10, 1778+10)

    # mat1.rebin("Ex", factor=5)
    # mat2.rebin("Ex", factor=5)
    # fig, ax = plt.subplots()
    # # mat1 /= mat2
    # values1, E = mat1.projection("Ex", Emin=1000, Emax=2000)
    # values2, E = mat2.projection("Ex", Emin=1000, Emax=2000)
    # ax.plot(E, values1/values2)
    # fig, ax = plt.subplots()

    # ax.plot(E, values1)
    # ax.plot(E, values2)
    # plt.show()

    return np.c_[E, values]


if __name__ == "__main__":
    # fwhm_pars = np.array([73.2087, 0.50824, 9.62481e-05])
    # Frank June 2020
    fwhm_pars = np.array([60.6499, 0.458252, 0.000265552])


    # # 144Nd -- Ex = 696 KeV
    # fname_exp = "exp/inbeam/144nd/alfna.m"
    # fname_bg = "exp/inbeam/144nd/alfna_bg.m"
    # Efit_low = 650
    # Efit_high = 720
    # Ecompare_low = 50
    # Ecompare_high = 1000
    # bg_ratio = 0

    # exp = arr_from_py(fname_exp, 640-40, 640+40, remove_negative=True)
    # bg = arr_from_py(fname_bg, 640-40, 640+40, remove_negative=True)
    # # fig, ax = plt.subplots()
    # # ax.plot(exp[:, 0], exp[:, 1], label="exp")
    # # # ax.plot(bg[:, 0], bg[:, 1], label="bg")

    # # # exp = arr_from_py(fname_exp, 640-20, 640+20)
    # # # ax.plot(exp[:, 0], exp[:, 1], label="exp2")

    # # # exp = arr_from_py(fname_exp, 640-10, 640+10)
    # # # ax.plot(exp[:, 0], exp[:, 1], label="exp3")

    # # ax.set_yscale("log")
    # # ax.legend()
    # # plt.show()

    # get_fom("696",
    #         exp, bg, fwhm_pars, bg_ratio,
    #         Efit_low, Efit_high,
    #         do_plot=True, printout=True)


    # 28Si -- Ex = 1.7 KeV
    fname_exp = "exp/inbeam/si-run1/28si/alfna.m"
    fname_bg = "exp/inbeam/si-run1/28si/alfna_bg.m"
    Efit_low = 1779-80
    Efit_high = 1779+35
    Ecompare_low = 50
    Ecompare_high = 1000
    bg_ratio = 0

    exp = arr_from_py(fname_exp, 1600-40, 1600+40, remove_negative=True)
    bg = arr_from_py(fname_bg, 1600-40, 1600+40, remove_negative=True)
    # exp = arr_from_py(fname_exp, 1516, 1575, remove_negative=True)
    # bg = arr_from_py(fname_bg,   1516, 1575, remove_negative=True)
    # exp = arr_from_py(fname_exp, 1575, 1625, remove_negative=True)
    # bg = arr_from_py(fname_bg,   1575, 1625, remove_negative=True)

    # small recalibration
    exp[:, 0] -= 2
    bg[:, 0] -= 2

    # test
    # bg[0, 1] /= 0

    # fig, ax = plt.subplots()
    # ax.plot(exp[:, 0], exp[:, 1], label="exp")
    # ax.plot(bg[:, 0], bg[:, 1], label="bg")

    # fig, ax = plt.subplots()
    # ax.plot(exp[:, 0], exp[:, 1], label="exp")
    # ax.plot(bg[:, 0], bg[:, 1], label="bg")

    # exp = arr_from_py(fname_exp, 1700-100, 1700+100, remove_negative=True)
    # ax.plot(exp[:, 0], exp[:, 1], label="exp2")

    # exp = arr_from_py(fname_exp, 1700-150, 1700+150, remove_negative=True)
    # ax.plot(exp[:, 0], exp[:, 1], label="exp3")

    # ax.set_yscale("log")
    # ax.legend()
    # plt.show()

    get_fom("1779",
            exp, bg, fwhm_pars, bg_ratio,
            Efit_low, Efit_high,
            do_plot=True, printout=True, xmax_plot=2200,
            ymin=1.2e2, ymax=1e5)

    """
    # 28Si -- Ex = ~4.6 KeV
    fname_exp = "exp/inbeam/si-run1/28si/alfna.m"
    fname_bg = "exp/inbeam/si-run1/28si/alfna_bg.m"
    Efit_low = 2750
    Efit_high = 2870
    Ecompare_low = 50
    Ecompare_high = 1000
    bg_ratio = 0

    exp = arr_from_py(fname_exp, 4300-40, 4300+40, remove_negative=True)
    bg = arr_from_py(fname_bg, 4300-40, 4300+40, remove_negative=True)

    # exp = arr_from_py(fname_exp, 3800, 4200, remove_negative=True)
    # bg = arr_from_py(fname_bg, 3800, 4200, remove_negative=True)

    # small recalibration
    exp[:, 0] -= 2
    bg[:, 0] -= 2

    fig, ax = plt.subplots()
    ax.plot(exp[:, 0], exp[:, 1], label="exp")
    ax.plot(bg[:, 0], bg[:, 1], label="bg")

    # exp = arr_from_py(fname_exp, 4300-100, 4300+100, remove_negative=True)
    # ax.plot(exp[:, 0], exp[:, 1], label="exp2")

    # exp = arr_from_py(fname_exp, 4300-150, 4300+150, remove_negative=True)
    # ax.plot(exp[:, 0], exp[:, 1], label="exp3")

    # ax.set_yscale("log")
    # ax.legend()

    from scipy.stats import norm
    from scipy.optimize import curve_fit

    def fitnorm(x, c, loc, scale, offset):
        return c * norm.pdf(x, loc, scale) + offset

    iE1 = np.abs(exp[:, 0] - 1720).argmin()
    iE2 = np.abs(exp[:, 0] - 1890).argmin()
    xfit = exp[iE1:iE2+1, 0]
    yfit = exp[iE1:iE2+1, 1]
    p0 = (20*np.max(yfit), 1778, 20, 400)
    bounds = (0, [np.inf, np.inf, np.inf, 410])

    popt, pcov = curve_fit(fitnorm, xfit, yfit, p0=p0,
                           bounds=bounds)
    ax.plot(xfit, fitnorm(xfit, *popt))

    iE1 = np.abs(exp[:, 0] - (2837-60)).argmin()
    iE2 = np.abs(exp[:, 0] - (2837+140)).argmin()
    xfit1 = exp[iE1:iE2+1, 0]
    yfit1 = exp[iE1:iE2+1, 1]
    p0 = (20*np.max(yfit1), 2837, 20, 5)

    bounds = (0, [np.inf, np.inf, np.inf, 10])

    popt1, pcov1 = curve_fit(fitnorm, xfit1, yfit1, p0=p0, bounds=bounds)
    ax.plot(xfit1, fitnorm(xfit1, *popt1))
    ax.set_yscale("log")

    fig, ax = plt.subplots()
    ax.plot(xfit, yfit-fitnorm(xfit, *popt))

    fig, ax = plt.subplots()
    ax.plot(xfit1, yfit1-fitnorm(xfit1, *popt1))

    plt.show()

    get_fom("4617",
            exp, bg, fwhm_pars, bg_ratio,
            Efit_low, Efit_high,
            do_plot=True, printout=True, xmax_plot=4800)
    """

    # # 12C -- Ex = 4.4 MeV
    fname_exp = "exp/inbeam/12C/h_bg_subtracted.m"
    # fname_bg = "exp/inbeam/si-run1/28si/alfna_bg.m"
    Efit_low = 4440-80
    Efit_high = 4440+35
    Ecompare_low = 50
    Ecompare_high = 1000
    bg_ratio = 0

    # especially 21, 26, 30 have "double peaks"
    # slight impresicion: at the moment not removed from simulations
    # but spectra there are almost the same for each detector

    good_dets = np.arange(30)

    # bad_dets = [3, 5, 7, 8, 18, 21, 26, 30]
    # bad_dets.sort()
    # bad_dets = np.array(bad_dets)
    # good_dets = \
    #     good_dets[bad_dets[np.searchsorted(bad_dets, good_dets)] != good_dets]

    exp_mat = Matrix(path=fname_exp)
    # exp_mat.plot()

    exp = np.zeros((len(exp_mat.Eg), 2))
    exp[:, 0] = exp_mat.Eg
    exp[:, 1] = exp_mat.values[good_dets, :].sum(axis=0)
    exp[:, 1][exp[:, 1] < 0] = 0  # remove negative counts for comparison

    exp_all = np.zeros((len(exp_mat.Eg), 2))
    exp_all[:, 0] = exp_mat.Eg
    exp_all[:, 1] = exp_mat.values[:, :].sum(axis=0)

    # dummy bg
    bg = exp.copy()

    fig, ax = plt.subplots()
    ax.plot(exp[:, 0], exp[:, 1], label="exp")
    ax.plot(exp_all[:, 0], exp_all[:, 1], label="all detectors")
    # ax.plot(bg[:, 0], bg[:, 1], label="bg")

    # fig, ax = plt.subplots()
    # ax.plot(exp[:, 0], exp[:, 1], label="exp")
    # ax.plot(bg[:, 0], bg[:, 1], label="bg")

    # exp = arr_from_py(fname_exp, 1700-100, 1700+100, remove_negative=True)
    # ax.plot(exp[:, 0], exp[:, 1], label="exp2")

    # exp = arr_from_py(fname_exp, 1700-150, 1700+150, remove_negative=True)
    # ax.plot(exp[:, 0], exp[:, 1], label="exp3")

    # ax.set_yscale("log")
    # ax.legend()
    # plt.show()

    get_fom("4440",
            exp, bg, fwhm_pars, bg_ratio,
            Efit_low, Efit_high,
            do_plot=True, printout=True, xmax_plot=5000)
