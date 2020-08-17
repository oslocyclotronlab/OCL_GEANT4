""" Run geant4 simulations on a grid for differernt geometry settings """
import numpy as np
from typing import Dict, List, Optional
from pathlib import Path


class MacroGen:
    def __init__(self):
        self.detsList = np.arange(30) + 1
        self.dgeometry = None

        self.run_macs = ["run_152Eu.mac", "run_60Co.mac", "run_133Ba.mac",
                         "run_137Cs.mac"]
        self.outdir = "../data/"
        self.outname_base = ""  # optionally: add eg sim_001_

        self._base_geometry_cmd = "/control/execute setup_rad_source.mac"
        pass

    def compose(self):
        return '\n'.join(*[self.geometry() + self.run()])

    def save(self, fname):
        fn = Path(fname)
        fn.write_text(self.compose())

    def geometry(self, dgeometry: Optional[Dict] = None,
                 unit: str = "cm") -> List[str]:
        d = self.dgeometry if dgeometry is None else dgeometry
        string = [
            self._base_geometry_cmd,
            f"/OCL/det/oscar/setCoatThickFront {d['front']:.4f} {unit}",
            f"/OCL/det/oscar/setCoatThickRad {d['radial']:.4f} {unit}",
            f"/OCL/det/oscar/setLaBrRefThick {d['reflector']:.4f} {unit}",
            f"/OCL/det/oscar/setLaBrLidHalfThick {d['lidhalf']:.4f} {unit}",
            "",
            # note: default unit here is cm!
            *[f"/OCL/det/oscar/setLaBrDist {N} {d['det']:.4f}"
              for N in self.detsList],
            "",
            "/run/initialize",
            ""]
        return string

    def run(self) -> List[str]:
        outname_base = Path(self.outdir)
        fnrun = [Path(run) for run in self.run_macs]
        fnout = [Path(self.outname_base + run) for run in self.run_macs]
        fnout = [outname_base / fn.with_suffix(".root") for fn in fnout]

        def basestring(fnrun, fnout):
            res = [f"/OCL/setOutName {fnout}",
                   f"/control/execute {fnrun}"]
            return res

        string = [basestring(fnrun_, fnout_)
                  for fnrun_, fnout_ in zip(fnrun, fnout)]
        flat_list = [item for sublist in string for item in sublist]
        return flat_list


def meshgrid_1array(*args):
    """ Make one grid of points from n different 1d grids.

    Example:
    [x1, x2] ,  [y1, y2, y3]
    -> [[x1, y1],
        [x2, y1],
        [x1, y2],
        [x2, y2],
        [x1, y3],
        [x2, y3]]
    """
    return np.array(np.meshgrid(*args)).T.reshape(-1, len(args))


if __name__ == "__main__":
    np.random.seed(65432)
    macgen = MacroGen()
    # print(macgen.run())

    # nominal thicknesses in cm
    dnominal = {'front': 0.2, 'radial': 0.1, 'reflector': 0.1,
                'lidhalf': 0.1, 'det': 16}
    dvar = {}
    for key, value in dnominal.items():
        dvar[key] = np.linspace(value*0.9, value*1.1, num=5)

    grid = meshgrid_1array(*dvar.values())
    print("Items to calculate: ", grid.shape)

    fnbase = Path("grid_macs")
    fnbase.mkdir(exist_ok=True)
    for i, pars in enumerate(grid):
        # print(f"Simulating gridpoint {i}")
        dtmp = {}
        for idict, key in enumerate(dnominal.keys()):
            dtmp[key] = pars[idict]

        macgen.outname_base = f"grid_{i}_"
        macgen.dgeometry = dtmp
        macro = macgen.save(fnbase / f"grid_{i}.mac")

    # create summary file with commands to run
    # due to concerns on calculation time we may not calculate all values,
    # but start with a random selection
    indices = np.arange(len(grid))
    np.random.shuffle(indices)

    cmds = [f"./OCL grid_macs/grid_{i}.mac" for i in indices]
    cmd_string = "\n".join(*[cmds])
    fn_sum = Path("grid_cmd_all.txt")
    fn_sum.write_text(cmd_string)




    # radial = np.linspace(1, 4, num=3)
    # z = np.linspace(10, 40, num=4)
    # xv = meshgrid_1array(x, y, z)
    # print(xv)
    # print(yv)

