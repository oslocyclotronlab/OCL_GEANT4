import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm

import ompy as om


def fFWHM(E, p):
    # sgima in keV
    return np.sqrt(p[0] + p[1] * E + p[2] * E**2)


def get_area(vec, E1, E2):
    return vec[vec.index(E1):vec.index(E2)+1].sum()


def get_counts_bgsubtract(vec, energy, extend=10, ext_bg=None):
    """ Get background subtracted nr counts at energy [+- extend]

    Args:
        vec: vector to process
        energy: photopeak energy
        extend: get area for energy+-extend, as
            counts are often placed in some bins around the "sharp" energy
        ext_bg: extend where bg should be taken into account

    Returns:
        peak_wo_bg, bg_avg: Peak with out background, and background
    """

    if ext_bg is None:
        ext_bg = extend+4

    bg_below = get_area(vec, energy-ext_bg, energy-extend)
    bg_above = get_area(vec, energy+extend,  energy+ext_bg)
    # note that bg is summed one additional channel
    bg_avg = (bg_below + bg_above) / (2*(ext_bg-extend+1))

    peak_counts = get_area(vec, energy-extend, energy+extend)
    peak_wo_bg = peak_counts - bg_avg*extend

    return peak_wo_bg, bg_avg


def get_efficiencies(vec, energy):
    """ Get total, photopeak (...) efficiencies

    Args:
        vec: vector to analyze, should not be folded with energy resolution yet
        energy: Energy of the FE peak

    """
    mama_vec = vec.copy()
    df = pd.Series()
    df["E"] = energy
    df["tot"] = vec.values.sum() / nevent
    df["tot>50keV"] = vec.values[vec.index(50):].sum() / nevent
    df["tot>200keV"] = vec.values[vec.index(200):].sum() / nevent
    df["tot>500keV"] = vec.values[vec.index(500):].sum() / nevent

    eff_photo, bg = get_counts_bgsubtract(vec, energy, extend=10) / nevent
    replace_counts(mama_vec, energy-5, energy+5, bg*nevent)
    df["fe"] = eff_photo

    if energy > 2*511:
        eff_, bg = get_counts_bgsubtract(vec, energy-511, extend=5) / nevent
        e = energy-511
        replace_counts(mama_vec, e-5, e+5, bg*nevent)
        df["se"] = eff_

        eff_, bg = get_counts_bgsubtract(vec, energy-511*2, extend=5) / nevent
        e = energy-511*2
        replace_counts(mama_vec, e-5, e+5, bg*nevent)
        df["de"] = eff_

        eff_, bg = get_counts_bgsubtract(vec, 511, extend=5) / nevent
        e = 511
        replace_counts(mama_vec, e-5, e+5, bg*nevent)
        df["511"] = eff_
    else:
        df["se"] = 0
        df["de"] = 0
        df["511"] = 0
    return df, mama_vec


def replace_counts(vec, E1, E2, counts):
    """ Replace counts in vec by new value """
    vec[vec.index(E1):vec.index(E2)+1] = counts


def efficiency_plots(efficiencies: pd.DataFrame, energy_grid):
    fig, ax = plt.subplots()
    # eff.plot(x="E", ax=ax)
    ax.plot(energy_grid, eff["tot>50keV"], "C0-", label="tot > 50 keV")
    ax.plot(energy_grid, eff["tot>200keV"], "C0--", label="tot > 200 keV")
    ax.plot(energy_grid, eff["tot>500keV"], "C0-.", label="tot > 500 keV")
    ax.plot(energy_grid, eff["fe"], "C1", label="full energy")
    ax.plot(energy_grid, eff["se"], "C1", label="single escape",
            ls=(0, (5, 5)))
    ax.plot(energy_grid, eff["de"], "C1", label="double escape",
            ls=(0, (3, 5, 1, 5)))
    ax.plot(energy_grid, eff["511"], "C2", label="annihilation",
            ls=(0, (3, 1, 1, 1)))
    ax.set_xlabel("Energy [keV]")
    ax.set_ylabel(r"efficiency $\epsilon$")
    ax.set_xlim(50, None)

    d = 16.3
    eff_geom = 30*(np.pi*(8.98/2)**2) / (4*np.pi*d**2)  # noqa
    ax.axhline(eff_geom, color="k", ls="--", alpha=0.5)

    ax.legend(loc="upper left", ncol=2)

    fig.tight_layout(pad=0.02)
    return fig, ax

def export_Guttormsen1996_interpolation():
    """ write files for Guttormsen1996 response interpolation

    TODO: clean up / refactor
    """
    for i, (energy, nevent) in enumerate(zip(tqdm(energy_grid), nevents)):
        fn = specdir / f"grid_{energy}keV_n{nevent}.root.m"
        if not fn.exists():
            continue

        vec = om.Vector(path=fn, units="MeV")
        vec.to_keV()

        eff.loc[i], mama_vec = get_efficiencies(vec, eff["E"][i])

        # rebin and smooth; rebin first to save time smoothing
        vec.rebin(mids=energy_out_uncut)
        vec.values = om.gauss_smoothing(vec.values, vec.E,
                                        fFWHM(vec.E, fwhm_pars))
        vec.rebin(mids=energy_out)

        respmat.values[i, :] = vec.values

        ompy_vec = mama_vec.copy()
        ompy_vec.cut(Emax=energy + 2000)
        ompy_vec.save(f"{fnompy}/cmp{energy}", filetype="mama", _assert=False)
        if energy < 10000:
            mama_vec.cut(Emax=10000)
            mama_vec.rebin(factor=2)
        else:
            mama_vec.rebin(factor=3)

        if energy != 50:
            mama_vec.save(f"{fnmama}/cmp{energy}", filetype="mama")

        eff_mama.loc[i] = eff.loc[i]
        eff_mama.loc[i, ["fe", "se", "de", "511"]] *= nevent
        # # get before smoothing, not after -> large effect
        # # eff_photo = get_area(vec, energy-10, energy+10) / nevent
        # # print(eff_photo)

        # vec.values /= i  # scale number of events
        # # vec.save(f"export/resp_m{i}.m")
        # vec.plot(ax=ax)

    eff.to_csv("efficiencies.csv")

    eff_mama = eff_mama.iloc[1:]
    eff_mama.drop(["tot>50keV", "tot>200keV", "tot>500keV"], axis=1,
                  inplace=True)
    fwhm_rel = fFWHM(eff_mama["E"], fwhm_pars) / eff_mama["E"]
    idx = abs(eff_mama['E'] - 1333).idxmin()
    fwhm_rel_1333 = fwhm_rel[idx]
    eff_mama.insert(1, "FWHM", fwhm_rel/fwhm_rel_1333)

    eff_mama["tot"] /= eff_mama["tot"][idx]
    mama_header =(
        "# OSCAR 2020 LaBr; 16.3 cm distance\n"
        "# Response functions from article arXiv:2008.06240v2 | Response functions simulated by GEANT4 (Fabio and Gry, August 2020) || FWHM rel. and Eff Normalized to 1 at 1.33 MeV.\n"
        "# Next: Numer of Lines\n"
        f"{len(eff_mama)}\n"
        "# Eg        FWHM_rel      Eff_tot       FE            SE            DE            c511"
        )
    for fdir in [fnmama, fnompy]:
        np.savetxt(fdir / "resp.dat", eff_mama.to_numpy(), fmt='%.4e',
                   delimiter="\t", header=mama_header, comments="")
        print("saved mama")


if __name__ == "__main__":
    figs_dir = Path("figs")
    figs_dir.mkdir(exist_ok=True)

    response_outdir = Path("response_export")
    response_outdir.mkdir(exist_ok=True)

    energy_grid = np.arange(50, 1e4, 10, dtype=int)
    nevents = np.linspace(6e5, 3e6, len(energy_grid), dtype=np.int)

    energy_grid = np.append(energy_grid, [int(1.2e4), int(1.5e4), int(2e4)])
    nevents = np.append(nevents, [int(3e6), int(3e6), int(3e6)])

    fwhm_pars = np.array([60.6499, 0.458252, 0.000265552])

    energy_out = np.arange(energy_grid[0], 21000, 10)
    energy_out_uncut = np.arange(0, 21000, 10)
    # response mat. with x: incident energy; y: outgoing
    respmat = om.Matrix(values=np.zeros((len(energy_grid), len(energy_out))),
                        Ex=energy_grid,
                        Eg=energy_out)

    eff = pd.DataFrame(columns=["E", "tot", "tot>50keV", "tot>200keV",
                                "tot>500keV",
                                "fe", "se", "de", "511"])
    eff["E"] = energy_grid

    eff_mama = eff.copy()

    specdir = Path("from_geant")
    # fig, ax = plt.subplots()
    # fig_plain, ax_plain = plt.subplots()

    fnmama = Path("mama_export")
    fnmama.mkdir(exist_ok=True)

    fnompy = Path("ompy_interpolate")
    fnompy.mkdir(exist_ok=True)

    # export files for interpolation with mama / ompy
    export_Guttormsen1996_interpolation()

    # different version of response matrix
    fn = "response_unnorm"
    _, ax, fig = respmat.plot(title="response unnormalized",
                              scale="log", vmin=1)
    ax.set_ylabel(r"$E_{incident}$ [KeV]")
    fig.savefig(figs_dir/f"{fn}.png")
    # respmat.save(response_outdir/f"{fn}.m") # too large for mama
    respmat.save(response_outdir/f"{fn}.txt")
    fig.clear()
    plt.close(fig)

    fn = "response_norm_efficiency"
    mat = respmat.copy()
    mat.values /= nevents[:, np.newaxis]
    _, ax, fig = mat.plot(title="response normalized to total eff.",
                          scale="log", vmin=1e-5, vmax=1e-1)
    ax.set_ylabel(r"$E_{incident}$ [KeV]")
    fig.savefig(figs_dir/f"{fn}.png")
    # mat.save(response_outdir/f"{fn}.m") # too large for mama
    mat.save(response_outdir/f"{fn}.txt")
    fig.clear()
    plt.close(fig)

    fn = "response_norm_1"
    mat = respmat.copy()
    mat.values /= mat.values.sum(axis=1)[:, np.newaxis]
    _, ax, fig = mat.plot(title="response normalized to 1",
                          scale="log", vmin=1e-5, vmax=1e-1)
    ax.set_ylabel(r"$E_{incident}$ [KeV]")
    fig.savefig(figs_dir/f"{fn}.png")
    # mat.save(response_outdir/f"{fn}.m") # too large for mama
    mat.save(response_outdir/f"{fn}.txt")
    fig.clear()
    plt.close(fig)

    fn = "response_squarecut_50keV_10.000keV_efficiency"
    mat = respmat.copy()
    mat.values /= nevents[:, np.newaxis]  # scale before energy cuts!
    mat.cut(axis="Ex", Emin=50, Emax=10e3)
    mat.cut(axis="Eg", Emin=50, Emax=10e3)
    _, ax, fig = mat.plot(title="response (square cut 50keV - 10 MeV\n"
                                "normalized to total efficiency",
                          scale="log", vmin=1e-5, vmax=1e-1)
    ax.set_ylabel(r"$E_{incident}$ [KeV]")
    fig.savefig(figs_dir/f"{fn}.png")
    mat.save(response_outdir/f"{fn}.m")
    mat.save(response_outdir/f"{fn}.txt")
    fig.clear()
    plt.close(fig)

    # comparison fig as spectra
    fig, ax = plt.subplots()
    for Ein in [0.5, 1, 2, 3, 5, 9]:
        mat.plot_projection(axis="Eg", ax=ax, Emin=Ein*1e3, Emax=Ein*1e3,
                            label=f"E_inc = {Ein} MeV",
                            scale="log")
    ax.set_ylim(1e-5, 1e-1)
    ax.set_ylabel("response/ bin /n_incident")
    ax.set_xlabel("Energy [keV]")
    ax.legend()
    fig.savefig(figs_dir/"response_functions.png")

    fn = "response_squarecut_50keV_10.000keV_norm_1"
    mat = respmat.copy()
    mat.cut(axis="Ex", Emin=50, Emax=10e3)
    mat.cut(axis="Eg", Emin=50, Emax=10e3)
    mat.values /= mat.values.sum(axis=1)[:, np.newaxis]
    _, ax, fig = mat.plot(title="response (square cut 50keV - 10 MeV\n"
                                "response normalized to 1",
                          scale="log", vmin=1e-5, vmax=1e-1)
    ax.set_ylabel(r"$E_{incident}$ [KeV]")
    fig.savefig(figs_dir/f"{fn}.png")
    mat.save(response_outdir/f"{fn}.m")
    mat.save(response_outdir/f"{fn}.txt")
    fig.clear()
    plt.close(fig)

    # efficiency plots
    fig, ax = efficiency_plots(eff, energy_grid)
    fig.savefig(figs_dir/"eff.png", dpi=200)

    # plt.show()
