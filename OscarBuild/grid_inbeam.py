import numpy as np
from pathlib import Path
import pandas as pd

from grid_simulations import MacroGen

if __name__ == "__main__":
    np.random.seed(65432)
    macgen = MacroGen()
    macgen._base_geometry_cmd = "/control/execute setup_normal_run.mac"
    macgen.run_macs = [#"run696keV.mac",
                       # "run1779keV.mac", "run2838keV.mac",
                       "runEx4617keV.mac",
                       "run4440keV.mac"
                       ]
    # print(macgen.run())

    # nominal thicknesses in cm
    dnominal = {'front': 0.2, 'radial': 0.1, 'reflector': 0.1,
                'lidhalf': 0.1, 'det': 16}
    dthick = {'front': 0.2, 'radial': 0.1, 'reflector': 0.1,
              'lidhalf': 0.2, 'det': 16}
    dsaintgb = {'front': 0.08, 'radial': 0.08, 'reflector': 0.22,
                'lidhalf': 0.1, 'det': 16}

    grid = [dnominal, dthick, dsaintgb]
    print("Items to calculate: ", len(grid))

    fnbase = Path("inbeam_grid_macs")
    fnbase.mkdir(exist_ok=True)
    for i, pars in enumerate(grid):
        # print(f"Simulating gridpoint {i}")
        dtmp = pars

        macgen.outname_base = f"inbeam_grid_{i}_"
        macgen.dgeometry = dtmp
        macro = macgen.save(fnbase / f"grid_{i}.mac")

    # create summary file with commands to run
    # due to concerns on calculation time we may not calculate all values,
    # but start with a random selection
    indices = np.arange(len(grid))
    # np.random.shuffle(indices)

    cmds = [f"./OCL inbeam_grid_macs/grid_{i}.mac" for i in indices]
    cmd_string = "\n".join(*[cmds])
    fn_sum = Path("inbeam_grid_cmd_all.txt")
    fn_sum.write_text(cmd_string)

    # grid_out = np.column_stack((np.arange(len(grid)), grid))
    # grid_out = pd.DataFrame(grid_out, columns=["grid_point", *dnominal.keys()])
    # grid_out = grid_out.astype({"grid_point": 'int'}, copy=False)
    # grid_out.to_pickle("inbeam_grid_inbeam.pickle")
