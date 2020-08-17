import numpy as np
from typing import Dict, List, Optional
from pathlib import Path

class MacroGenResponse:
    def __init__(self, energy: Optional[float] = None,
                 nevent: Optional[int] = None):
        self.energy = energy
        self.nevent = nevent
        self.outdir = "../data/"
        self.outname_base = "grid_"  # optionally: add eg sim_001_

        self._base_geometry_cmd = "/control/execute setup_normal_run.mac"

    def compose(self):
        return '\n'.join(*[self.geometry() + self.run()])

    def save(self, fname):
        fn = Path(fname)
        fn.write_text(self.compose())

    def geometry(self,
                 unit: str = "cm") -> List[str]:
        string = [
            self._base_geometry_cmd,
            "/run/initialize",
            ""]
        return string

    def run(self) -> List[str]:
        assert np.issubdtype(type(self.nevent), np.integer)

        outname_base = Path(self.outdir)
        fnout = Path(f"{self.outname_base}{self.energy}keV_n{self.nevent}")
        fnout = outname_base / fnout.with_suffix(".root")

        def basestring(energykeV, fnout, nevent):
            res = ["# Particle type, position, energy...",
                   "/gps/particle gamma",
                   "/gps/number 1",
                   "",
                   "# Particle source distribution",
                   "/gps/pos/type Plane",
                   "/gps/pos/shape Ellipse",
                   "/gps/pos/centre 0. 0. 0. mm",
                   "/gps/pos/halfx 0.75 mm",
                   "/gps/pos/halfy 1.25 mm",
                   "/gps/ang/type iso",
                   "",
                   f"/gps/energy {energykeV} keV",
                   f"/OCL/setOutName {fnout}",
                   "",
                   "# Number of events to run",
                   f"/run/beamOn {nevent}",
                   ]
            return res

        string = basestring(self.energy, fnout, self.nevent)
        # flat_list = [item for sublist in string for item in sublist]
        return string


if __name__ == "__main__":
    energy_grid = np.arange(50, 1e4, 10, dtype=int)
    nevents = np.linspace(6e5, 3e6, len(energy_grid), dtype=np.int)

    energy_grid = np.append(energy_grid, [int(1.2e4), int(1.5e4), int(2e4)])
    nevents = np.append(nevents, [int(3e6), int(3e6), int(3e6)])

    fnbase = Path("response_grid_macs")
    fnbase.mkdir(exist_ok=True)
    for i, (energy, nevent) in enumerate(zip(energy_grid, nevents)):
        # print(f"Simulating gridpoint {i}")
        macgen = MacroGenResponse(energy=energy, nevent=nevent)
        macro = macgen.save(fnbase / f"grid_{i}.mac")

    # create summary file with commands to run
    # sorted by decreasing computation time (highest energies first)
    indices = np.arange(len(energy_grid))
    # np.random.shuffle(indices)

    cmds = [f"./OCL {fnbase}/grid_{i}.mac > $LOGDIR/out.o$LAUNCHER_JID"
            for i in indices[::-1]]
    cmd_string = "\n".join(*[cmds])
    fn_sum = Path("response_grid_cmds.txt")
    fn_sum.write_text(cmd_string)
