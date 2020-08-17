# Response functions generated from the Oscar OCL model

## General structure
- It is assumed that the response functions are processed from the root trees
  to single histogram for all detectors combined, see `from_geant` folder. One
  might for example use `export_hist.C` for that. The naming convention is
  `grid_{incidentEnergy}keV_n{numberEventsSimulated}.root.m`.
  Note that the files might be too large (too many bins) for mama to read propperly,
  the postprocessing was done with OMpy.
- Postprocessing in `spec_to_matrix.py`: Gets efficiencies, folds response functions
  and generates response matrices
- `figs/`: generated figures
- `response_export`: Exported response functions, several formats
- `efficiencies.csv`: efficiency results from `spec_to_matrix.py`
- The `compare[...]` scripts were used to compare the experimental
  spectra to the geant4 simulations. The code is not *cleaned up*.


## Naming convention for response matrices:
- `_unnorm`: Plain response function as calculated from geant4 
  (however, here and below, folded with energy response) 
- `_norm_efficiency`: spectra for each incident energy are normalized to the 
  total efficiency (equals division by number of incident events)
- `_norm_1`: spectra for each incident energy are normalized to 1
- `squarecut_50keV_10.000keV`: Square response matrix, with a lower (higher) cut 
   of 50 keV and 10 MeV, respectively. These are also small enough to be 
   able to load them into mama.
