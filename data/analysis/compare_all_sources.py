import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
import fnmatch
import os
import re
from tqdm import tqdm
from datetime import datetime
import pandas as pd

from uncertainties import unumpy, ufloat, correlated_values
import sys

import pandas as pd
# sys.path.append("../exp")

from compare import SpectrumComparison, fFWHM, PeakFitter

from scipy.optimize import curve_fit

def save_coords_from_click(fig, fname="coords.txt"):
    try:
        coords = np.loadtxt(fname)
        print("loading file (instead)")
        # plt.show()
    except OSError:
        coords = []

        def onclick(event):
            """ depends on global `coords` """
            ix, iy = event.xdata, event.ydata
            print(f'x = {ix:.0f}, y = {iy:.0f}')
            coords.append((ix, iy))
            return coords

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()
        fig.canvas.mpl_disconnect(cid)

        np.savetxt(fname, np.array(coords),
                   header="x y")
    coords = np.array(coords)
    coords = coords.astype(int)
    return coords



def get_fom(fnisotope,
            fname_exp, fname_bg, fwhm_pars,
            measure_time_exp, measure_time_bg, idets,
            Efit_low, Efit_high,
            do_plot=True, printout=False,
            manual_ratio=None, fitpeaks=None,
            xmax=1600, figtext=None):
    """ get figure of merrit

    fnisotope: str, like "60Co", or "152Eu"
    """
    fname_sims = []
    for file in os.listdir('mama_spectra/root_files'):
        if fnmatch.fnmatch(file, f'grid_-5_*{fnisotope}*_all.m'):
            print("machting:", file)
            fname_sims.append(os.path.join("mama_spectra/root_files", file))
    fname_sims.sort()
    grid_points = np.full_like(fname_sims, np.nan, dtype=float)
    foms = np.zeros((len(fname_sims), 3))

    for i, fname_sim in enumerate(tqdm(fname_sims)):
        if i > 6:
            break
        if printout:
            print("fitting: ", fname_sim)

        grid_point = int(re.search(r"grid_(-\d*)_", fname_sim).groups()[0])
        grid_points[i] = grid_point

        sc = SpectrumComparison()
        sc.get_data(fname_sim, fname_exp, fname_bg, fwhm_pars,
                    measure_time_exp, measure_time_bg, idet=idets,
                    recalibrate=True)

        if manual_ratio is None:
            sc.scale_sim_to_exp_area(Efit_low, Efit_high)
        else:
            sc.scale_sim_to_exp_manual(manual_ratio)

        fit_peaks(fitpeaks, fnisotope, sc, grid_point)

        # compare_plain_counts(sc, manual_ratio, Elimits)
        # difference_total(sc)

        sc.get_chi2()
        rel_diff, rel_diff_smooth = sc.get_rel_diff(smooth_window_keV=20)

        foms[i, :] = sc.fom(Ecompare_low, Ecompare_high, printout=False)

        if do_plot:
            fig, (ax1, ax2) = sc.plots(title=fname_sim, xmax=xmax,
                                       plot_smoothed=False)

            ax1.text(0.5, 0.9, figtext, horizontalalignment='center',
                     verticalalignment='center', transform=ax1.transAxes,
                     fontsize="large")

            fig.savefig(f"figs/{fnisotope}_{grid_point:.0f}_noleg.png")

            ax1.legend(loc="best")
            # ax2.legend(loc="best")

            fig.savefig(f"figs/{fnisotope}_{grid_point:.0f}.png")

            add_peaks_plot(ax1, fnisotope, fitpeaks, grid_point)

            ax1.legend()
            fig.savefig(f"figs/{fnisotope}_{grid_point:.0f}_fits.png")
            # plt.show()
            plt.close(fig)

    if printout:
        ltab = [[name, *foms[i, :]] for i, name in enumerate(fname_sims)]
        print("\nComparisons between {} and {} keV:"
              .format(Ecompare_low, Ecompare_high))
        print(tabulate(ltab,
                       headers=["Name", "chi2", "rel_diff[%]",
                                "rel_diff_smoothed[%]"],
                       floatfmt=".2f"))
    df = pd.DataFrame(foms, columns=[f"chi2_{fnisotope}",
                                     f"rel_diff_{fnisotope}",
                                     f"rel_diff_smoothed_{fnisotope}"])
    df["grid_point"] = grid_points
    df = df[df.grid_point.notnull()]  # workaround if going through whole loop
    df = df.astype({"grid_point": 'int'}, copy=False)

    return df


def fit_peaks(fitpeaks, fnisotope, sc, grid_point):
    """ Fit peaks by different functions

    Note: changes dictrionary fitpeaks inplace!
    """
    fitpeak = fitpeaks[fnisotope] if fitpeaks is not None else {}
    for peak, specs in fitpeak.items():
        pf = PeakFitter(sc.exp[:, 0], sc.exp[:, 1],
                        unumpy.std_devs(sc.uyexp))
        _fit_peaks(pf, specs, f"exp_{grid_point}", grid_point)

        pf = PeakFitter(sc.sim[:, 0],
                        unumpy.nominal_values(sc.uysim_scaled),
                        unumpy.std_devs(sc.uysim_scaled))
        _fit_peaks(pf, specs, f"sim_{grid_point}", grid_point)


def _fit_peaks(pf, specs, key, grid_point):
    specs[key] = {}
    if isinstance(specs["E"], float):
        popt, pcov = pf.fitGausBgStep(specs["pre_region"],
                                      specs["post_region"])
        print(popt)
        print(pcov[0, 0])

        specs[key]["area"] = ufloat(popt[0], np.sqrt(pcov[0, 0]))

    elif len(specs["E"]) == 2:
        popt, pcov = pf.fitDoubleGausBgStep(specs["pre_region"],
                                            specs["post_region"],
                                            p0_E=specs["E"])
        specs[key]["area"] = [ufloat(popt[0], np.sqrt(pcov[0, 0])),
                              ufloat(popt[3], np.sqrt(pcov[3, 3]))]
    specs[key]["popt"] = popt
    specs[key]["xfit"] = pf._xfit


def add_peaks_plot(ax, fnisotope, fitpeaks, grid_point):
    fitpeak = fitpeaks[fnisotope] if fitpeaks is not None else {}
    for peak, specs in fitpeak.items():
        for item in [f"exp_{grid_point}", f"sim_{grid_point}"]:
            x = specs[item]["xfit"]
            popt = specs[item]["popt"]
            if isinstance(specs["E"], float):
                ffit = PeakFitter.gaus_bg_step
            elif len(specs["E"]) == 2:
                ffit = PeakFitter.doublegaus_bg_step
            ax.plot(x,
                    ffit(x, *popt), "--", label=f"fit_{peak}_{item}")


def compare_plain_counts(sc, manual_ratio, Elimits):
    """ compares counts within Elimits directly """
    ncounts_Elims = np.zeros((len(Elimits), 2))
    for j, (E1, E2) in enumerate(Elimits):
        ncounts_Elims[j, 0] = sc.get_area(sc.exp, E1, E2)
        sim_scaled = np.c_[sc.xsim, unumpy.nominal_values(sc.uysim_scaled)]
        ncounts_Elims[j, 1] = sc.get_area(sim_scaled, E1, E2)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    fig.suptitle(f"grid_point {grid_point}")
    counts_exp = unumpy.uarray(ncounts_Elims[:, 0],
                               np.sqrt(ncounts_Elims[:, 0]))
    counts_sim = unumpy.uarray(ncounts_Elims[:, 1]/manual_ratio,
                               np.sqrt(ncounts_Elims[:, 1])/manual_ratio)
    counts_sim *= manual_ratio
    counts_sim *= 5  # *5 due to binwidth
    diff = (counts_exp-counts_sim)/counts_exp * 100
    ax1.errorbar(Elimits.mean(axis=1), unumpy.nominal_values(counts_exp),
                 yerr=unumpy.std_devs(counts_exp), fmt="o")
    ax1.errorbar(Elimits.mean(axis=1), unumpy.nominal_values(counts_sim),
                 yerr=unumpy.std_devs(counts_sim), fmt="x")
    ax2.errorbar(Elimits.mean(axis=1), unumpy.nominal_values(diff),
                 yerr=unumpy.std_devs(diff), fmt="x")

    ax2.axhline(0, c="k", ls="--")
    ax2.set_ylim(-20, 20)
    return fig, (ax1, ax2)


def difference_total(sc):
    """ total differences between sim and exp in some area """
    sim_scaled = np.c_[sc.xsim, unumpy.nominal_values(sc.uysim_scaled)]
    tot_exp = sc.get_area(sc.exp, 100, 1300)
    tot_sim = sc.get_area(sim_scaled, 100, 1300)
    tot_sim *= 5  # 5 due to binwitdth
    tot_diff = (tot_exp-tot_sim)/tot_exp * 100
    print(f"Tot diff Elims1 [%]: {tot_diff:.2f}")

    tot_exp = sc.get_area(sc.exp, 50, 200)
    tot_sim = sc.get_area(sim_scaled, 50, 200)
    tot_sim *= 5  # 5 due to binwitdth
    tot_diff = (tot_exp-tot_sim)/tot_exp * 100
    print(f"Tot diff Elims2 [%]: {tot_diff:.2f}")


def ph_efficiency(x, p):
    """ also from gf3 """
    lneff_low = ln_pheff_per_region(x, p[0:3], Ec=100)
    lneff_high = ln_pheff_per_region(x, p[3:6], Ec=1000)
    inner = lneff_low**(-p[6]) + lneff_high**(-p[6])
    return np.exp(inner**(-1/p[6]))


def ph_efficiency_low(x, *p):
    """ also from gf3 """
    lneff_low = ln_pheff_per_region(x, *p, Ec=1)
    return np.exp(lneff_low)

def ph_efficiency_low_unumpy(x, *p):
    """ also from gf3 """
    lneff_low = ln_pheff_per_region(x, *p, Ec=1)
    return unumpy.exp(lneff_low)


def ln_pheff_per_region(x, p0, p1, p2, Ec=1):
    """ p0 + p1 * (x/Ec) + p2 * (x/Ec)**2 on log log plot"""
    logeff = p0 + p1*np.log(x/Ec) + p2*np.log(x/Ec)**2
    return logeff

def plot_photoeff(fitpeaks, grid_point, xmax=1500):
    fmt_list = ["<", ">", "^", "o", "x", "+", "."]
    fig, ax = plt.subplots(dpi=200)
    fig.suptitle(f"grid_point {grid_point}")

    effs_exp = []
    effs_sim = []
    # for name, group in grouped:
    #     ax.errorbar(group["E"]+3, group["eff"], yerr=group["sigma_eff"],
    #                 label=f"{name}_Frank_jitter", fmt="o", mfc="None",
    #                 alpha=0.2)

    for i_isotope, (isotope, disotope) in enumerate(fitpeaks.items()):
        try:
            ndecays = files[isotope]["ndecays"]
        except KeyError:
            continue
        if i_isotope < len(fmt_list):
            fmt = fmt_list[i_isotope]
        else:
            fmt = fmt_list[len(fmt_list)%i_isotope]

        for i, (peak, values) in enumerate(disotope.items()):
            try:
                print("before")
                val = values[f"exp_{grid_point}"]
                print("after")

                print(peak, values["E"])
                label = f"{isotope}" if i == 0 else None
                if isinstance(values["E"], float):
                    eff = val["area"]/values["intensity"] / ndecays / frac_dets

                    disotope[peak]["eff_exp"] = eff
                    ax.errorbar(values["E"], eff.nominal_value,
                                yerr=eff.std_dev,
                                c="C0", label=label, fmt=fmt, mfc="None")

                elif len(values["E"]) == 2:
                    eff = [area/intensity / ndecays / frac_dets
                           for area, intensity in zip(val["area"],
                                                      values["intensity"])]

                    disotope[peak]["eff_exp"] = eff
                    ax.errorbar(values["E"],
                                [p.nominal_value for p in eff],
                                yerr=[p.std_dev for p in eff],
                                c="C0", label=label, fmt=fmt, mfc="None")

                val = values[f"sim_{grid_point}"]

                # label = f"sim_{isotope}" if i == 0 else None
                label=None
                if isinstance(values["E"], float):
                    eff = (val["area"]/values["intensity"] / ndecays
                           / frac_dets)

                    disotope[peak]["eff_sim"] = eff
                    ax.errorbar(values["E"], eff.nominal_value,
                                yerr=eff.std_dev,
                                c="C1", label=label, fmt=fmt, mfc="None")
                elif len(values["E"]) == 2:
                    eff = [area/intensity / ndecays / frac_dets
                           for area, intensity in zip(val["area"],
                                                      values["intensity"])]

                    disotope[peak]["sim_exp"] = eff
                    ax.errorbar(values["E"],
                                [p.nominal_value for p in eff],
                                yerr=[p.std_dev for p in eff],
                                c="C1", label=label, fmt=fmt,
                                mfc="None")
                effs_exp.append([values["E"], disotope[peak]["eff_exp"]])
                effs_sim.append([values["E"], disotope[peak]["eff_sim"]])
            except KeyError:
                print("passing this:", grid_point, peak, values.keys())

    arr = np.array([[energy, y.nominal_value, y.std_dev]
                    for (energy, y) in effs_exp])
    arr = arr[arr[:, 0].argsort()]
    x = arr[:, 0]
    p0 = [60, 0.5, 2e-4]
    popt, pcov = curve_fit(ph_efficiency_low, x, arr[:, 1], p0=p0,
                           sigma=arr[:, 2])
    print(f"popt exp gp{grid_point}:", popt)
    print(f"pcov exp gp{grid_point}:\n",
          tabulate(np.c_[[r"$p_0$", r"$p_1$", r"$p_2$"], pcov],
                   tablefmt="latex_booktabs", floatfmt=".1e",
                   headers=[r"$p_0$", r"$p_1$", r"$p_2$"]))
    x = np.linspace(200, xmax, num=100)

    poptcorr = correlated_values(popt, pcov)
    print(f"popt exp gp: {poptcorr}")
    y = ph_efficiency_low_unumpy(x, *poptcorr)
    nom = unumpy.nominal_values(y)
    std = unumpy.std_devs(y)
    ax.plot(x, nom, "C0-", alpha=0.7)
    ax.fill_between(x, nom-std, nom+std, color="C0", alpha=0.15,
                    label=r"fit to $\epsilon_\mathrm{fe, exp}$")

    arr = np.array([[energy, y.nominal_value, y.std_dev]
                    for (energy, y) in effs_sim])
    arr = arr[arr[:, 0].argsort()]
    x = arr[:, 0]
    p0 = [60, 0.5, 2e-4]
    popt, pcov = curve_fit(ph_efficiency_low, x, arr[:, 1], p0=p0,
                           sigma=arr[:, 2])
    print(f"popt sim gp{grid_point}:", popt)
    print(f"pcov sim gp{grid_point}:\n",
          tabulate(np.c_[[r"$p_0$", r"$p_1$", r"$p_2$"], pcov],
                   tablefmt="latex_booktabs", floatfmt=".1e",
                   headers=[r"$p_0$", r"$p_1$", r"$p_2$"]))
    x = np.linspace(200, xmax, num=100)

    poptcorr = correlated_values(popt, pcov)
    print(f"popt sim gp: {poptcorr}")
    y = ph_efficiency_low_unumpy(x, *poptcorr)
    nom = unumpy.nominal_values(y)
    std = unumpy.std_devs(y)
    ax.plot(x, nom, "C1-", alpha=0.7)
    ax.fill_between(x, nom-std, nom+std, color="C1", alpha=0.15,
                    label=r"fit to $\epsilon_\mathrm{fe, sim}$")

    ax.legend()
    ax.set_xlabel("Energy [keV]")
    ax.set_ylabel(r"Full-energy peak efficiency $\epsilon_\mathrm{fe}$")
    return fig, ax

if __name__ == "__main__":
    # fwhm_pars = np.array([73.2087, 0.50824, 9.62481e-05])
    # Frank June 2020
    fwhm_pars = np.array([60.6499, 0.458252, 0.000265552])


    # print(fFWHM(80, fwhm_pars))
    # sys.exit()

    files = {
        "133Ba": {"t": 1049.48, "ndecays": 4.28917e+07},
        "60Co": {"t":  1123.57, "ndecays": 2.74162e+07},
        "152Eu": {"t": 1065.10, "ndecays": 3.45591e+07},
        "137Cs": {"t": 676.307, "ndecays": 2.77619e+07},
        "241Am": {"t": 969.774},
        "Bg": {"t":    1432.19, }}

    for key, file in files.items():
        try:
            file["ndecays"] = ufloat(file["ndecays"], 0.03*file["ndecays"])
        except KeyError:
            continue

    ndecays_sim = 2e5
    diff_binwidth = 5

    # dets with better resolution (internat peak structure vissible)
    idets = [1, 2, 6, 8, 10, 11, 12,
             14, 15, 16, 17, 18, 19, 20, 21, 22,
             24, 25, 27, 29]
    frac_dets = len(idets)/30

    # # 60Co
    Elimits = np.array([
                    [1110,  1230],
                    [1270,  1350]  # very short due to Bg subraction problem
                    ])

    fitpeaks = {name: {} for name in files.keys()}
    peak_energy = [1173.228, 1332.492]
    peak_intensity = [.9985, .999826]
    peak = {"E": peak_energy[0], "intensity": peak_intensity[0],
            "pre_region": [1075, 1100], "post_region": [1240, 1270]}
    fitpeaks["60Co"][int(peak["E"])] = peak
    peak = {"E": peak_energy[1], "intensity": peak_intensity[1],
            "pre_region": [1235, 1270], "post_region": [1400, 1470]}
    fitpeaks["60Co"][int(peak["E"])] = peak

    fname_exp = "exp/60Co.txt"
    fname_bg = "exp/Bg.txt"
    measure_time_exp = files["60Co"]["t"]  # //seconds
    measure_time_bg = files["Bg"]["t"]  # //seconds
    Efit_low = 1173 - 50 - 50
    Efit_high = 1173 + 50 + 50
    Ecompare_low = 50
    Ecompare_high = 2000

    df = get_fom("60Co",
                 fname_exp, fname_bg, fwhm_pars,
                 measure_time_exp, measure_time_bg, idets,
                 Efit_low, Efit_high,
                 do_plot=True, printout=True,
                 manual_ratio=files["60Co"]["ndecays"].nominal_value/ndecays_sim/diff_binwidth,
                 fitpeaks=fitpeaks,
                 xmax=1500,
                 figtext=r"$^{60}$Co")

    # df_all = df

    # 152Eu
    fname_exp = "exp/152Eu.txt"
    fname_bg = "exp/Bg.txt"
    measure_time_exp = files["152Eu"]["t"]  # //seconds
    measure_time_bg = files["Bg"]["t"]  # //seconds
    Efit_low = 720
    Efit_high = 830
    Ecompare_low = 50
    Ecompare_high = 1000

    # for "direct" comparison of counts (no peakfitting...)
    Elimits = np.array([
                    [105,  130],
                    [220,  260],
                    [315,  360],
                    [390,  420],
                    [420,  460],
                    [660,  700],
                    [740,  810],
                    [830,  900],
                    [910,  1000],
                    [1030, 1160],
                    [1170, 1240],
                    [1260, 1342]])

    # Photo-efficiency fits
    # Energy, intensity, pre_Elow, pre_Ehigh, post_Elow, post_Ehigh
    peaks = np.array(
        [[244.6974,  0.0755, 200 , 220 , 260 ,  280 ],   # noqa
         [344.2785,  0.2659, 300 , 315 , 378 ,  391 ],   # noqa
         [778.9045,  0.1293, 720 , 740 , 805 ,  830 ],   # noqa
         [867.38  ,  0.0423, 810 , 830 , 900 ,  930 ],   # noqa
         [964.057 ,  0.1451, 900 , 920 , 995 ,  1005],   # noqa
         [1408.013,  0.2087, 1340, 1350, 1550,  1650]  # bad fit
         ])  # noqa

    def peak_dict_from_arr(arr):
        peak = {"E": arr[0], "intensity": arr[1],
                "pre_region": [arr[2], arr[3]],
                "post_region": [arr[4], arr[5]]}
        return peak

    for i in range(len(peaks)):
        peak = peak_dict_from_arr(peaks[i, :])
        fitpeaks["152Eu"][int(peak["E"])] = peak

    def doublepeak_dict_from_arr(arr):
        peak = {"E": arr[:, 0], "intensity": arr[:, 1],
                "pre_region": [arr[0, 2], arr[0, 3]],
                "post_region": [arr[0, 4], arr[0, 5]]}
        return peak

    # # Energy, intensity, pre_Elow, pre_Ehigh, post_Elow, post_Ehigh
    # dobulepeaks = np.array(
    #     [[411.1165,  .02237, 380 , 390 , 460 ,  470 ],   # noqa
    #      [443.9606,  .02827, np.nan, np.nan, np.nan, np.nan],   # noqa
    #     ])  # noqa
    # fitpeaks["152Eu"]["411_444_double"] = doublepeak_dict_from_arr(dobulepeaks)

    df = get_fom("152Eu",
                 fname_exp, fname_bg, fwhm_pars,
                 measure_time_exp, measure_time_bg, idets,
                 Efit_low, Efit_high,
                 do_plot=True, printout=True,
                 manual_ratio=files["152Eu"]["ndecays"].nominal_value/ndecays_sim/diff_binwidth,
                 fitpeaks=fitpeaks,
                 xmax=1600,
                 figtext=r"$^{152}$Eu")
    # # df_all = df_all.merge(df, on="grid_point", how="outer")

    # 133Ba
    fname_exp = "exp/133Ba.txt"
    fname_bg = "exp/Bg.txt"
    measure_time_exp = files["133Ba"]["t"]  # //seconds
    measure_time_bg = files["Bg"]["t"]  # //seconds
    Efit_low = 285
    Efit_high = 320
    # Efit_low = 330
    # Efit_high = 400
    Ecompare_low = 50
    Ecompare_high = 300

    # # Energy, intensity, pre_Elow, pre_Ehigh, post_Elow, post_Ehigh
    # dobulepeaks = np.array(
    #     [[276.3989,  .0716, 245 , 255 , 320 ,  330],   # noqa
    #      [302.8508,  .1834, np.nan, np.nan, np.nan, np.nan],   # noqa
    #     ])  # noqa
    # fitpeaks["133Ba"]["276_302_double"] = doublepeak_dict_from_arr(dobulepeaks)

    # # Energy, intensity, pre_Elow, pre_Ehigh, post_Elow, post_Ehigh
    # dobulepeaks = np.array(
    #     [[356.0129, .6205, 320 , 330 , 410 ,  425],   # noqa
    #      [383.8485, .0894, np.nan, np.nan, np.nan, np.nan],   # noqa
    #     ])  # noqa
    # fitpeaks["133Ba"]["276_302_double"] = doublepeak_dict_from_arr(dobulepeaks)


    df = get_fom("133Ba",
                 fname_exp, fname_bg, fwhm_pars,
                 measure_time_exp, measure_time_bg, idets,
                 Efit_low, Efit_high,
                 do_plot=True, printout=False,
                 manual_ratio=files["133Ba"]["ndecays"]/ndecays_sim/diff_binwidth,
                 xmax=460,
                 figtext=r"$^{133}$Ba"
                 #fitpeaks=fitpeaks
                 )
    # # df_all = df_all.merge(df, on="grid_point", how="outer")

    # # 137Cs
    fname_exp = "exp/137Cs.txt"
    fname_bg = "exp/Bg.txt"
    measure_time_exp = files["137Cs"]["t"]  # //seconds
    measure_time_bg = files["Bg"]["t"]  # //seconds
    Efit_low = 600
    Efit_high = 700
    Ecompare_low = 50
    Ecompare_high = 300

    peak = {"E": 661.657, "intensity": .851,
            "pre_region": [580, 615], "post_region": [705, 740]}
    fitpeaks["137Cs"][int(peak["E"])] = peak

    df = get_fom("137Cs",
                 fname_exp, fname_bg, fwhm_pars,
                 measure_time_exp, measure_time_bg, idets,
                 Efit_low, Efit_high,
                 do_plot=True, printout=True,
                 manual_ratio=files["137Cs"]["ndecays"].nominal_value/ndecays_sim/diff_binwidth,
                 fitpeaks=fitpeaks,
                 xmax=800,
                 figtext=r"$^{137}$Cs")

    df = pd.read_csv('exp/labr_eff_fit_export.dat')
    new_name = df.columns[0].split("# ")[1]
    df = df.rename(columns={df.columns[0]: new_name})
    grouped = df.groupby("Isotope")

    for grid_point in [-5]:
        fig, ax = plot_photoeff(fitpeaks, grid_point)
        fig.savefig(f"figs/photoeff_{grid_point}.png")

    for grid_point in [-5]:
        fig, ax = plot_photoeff(fitpeaks, grid_point, xmax=15e3)

        # photopeak eff from geant (direct, not fitted)
        df = pd.read_csv("response/efficiencies.csv")
        df.plot(ax=ax, x="E", y="fe", label="geant4", color="k",
                linestyle="--")

        # Fits from Wanja, 6 Aug 2020'
        # 12C
        eff = {"E": [1778.969, 2838.29, 3200.7],
               "eff": [0.13462, 0.11538, 0.1053],
               "dE": [0.011, 0.15, 0.5],
               "deff": [0.13462*0.135, 0.11538*0.135, 0.1053*0.135]}
        ax.errorbar(x=eff["E"], xerr=eff["dE"], y=eff["eff"], yerr=eff["deff"],
                    label="28Si, Wanja", fmt="C3o")
        ax.errorbar(x=4439.82, xerr=0.21, y=0.08262, yerr=0.08262*0.135,
                    label="12C, Wanja", fmt="C4o")

        # Franks fits
        popt = [-0.29138, -0.0283581, -0.0271264]
        pcov = [[0.0135246, -0.00290798, 0.000121002, ],
                [-0.00290798, 0.000922994, -7.2436e-05, ],
                [0.000121002, -7.2436e-05, 8.42127e-06, ]]
        poptcorr = correlated_values(popt, pcov)
        x = np.linspace(200, 15e3, num=100)
        y = ph_efficiency_low_unumpy(x, *poptcorr)
        nom = unumpy.nominal_values(y)
        std = unumpy.std_devs(y)
        ax.plot(x, nom, "C4-", alpha=0.7)
        ax.fill_between(x, nom-std, nom+std, color="C4", alpha=0.15,
                        label=r"Frank's fit to $\epsilon_\mathrm{fe, exp}$")

        ax.legend()
        fig.savefig(f"figs/photoeff_{grid_point}_wtih_geant.png")
        ax.set_xlim(None, 5e3)
        ax.set_ylim(None, 0.32)
        fig.savefig(f"figs/photoeff_{grid_point}_wtih_geant_zoom.png")

    plt.show()
    # # df_all = df_all.merge(df, on="grid_point", how="outer")

    # now = datetime.now()
    # df_all.to_pickle(f'chi2_df_{now.strftime("%y%d%m_%H%M%S")}.pickle')
    # print(df_all[:8])
