# Atmospheric Radiative Transfer Modeling

This repository contains code and analysis for modeling radiative transfer phenomena in an idealized atmosphere model. The computational project models effects like refractive bending, scattering, and intensity changes as a function of wavelength and solar geometry.

## Overview

The Python code in `cp2.py` implements analytical radiative transfer equations for a multi-layer atmosphere model. Key phenomena modeled include:

- Refractive bending of ray path through layers
- Extinction due to Rayleigh scattering
- Scattering source terms 
- Reflection at the surface

The code computes the upwelling and downwelling radiation fractions over a range of solar zenith angles from 0 to 80 degrees. Effects on two wavelengths - 350nm and 1000nm - are analyzed.

## Contents

- `cp2.py`: Main Python code for the radiative transfer model
- `atmosphere_model.txt`: Layered atmosphere profile used for the model
- `report.pdf`: Writeup describing the modeling approach, results, and analysis
- `ourput/`: Directory containing plots generated from the code

## Usage

The code requires NumPy, SciPy, and Matplotlib. To run the simulation:

```
cp2.py
```

The model atmosphere and solar angles can be configured in the script.

## Results

Key results include:

- Refraction causes minimal path length change, more noticeable at higher solar zenith angles
- Upwelling fraction decreases exponentially with increasing solar zenith due to increased scattering
- Shorter wavelengths have lower downwelling radiation deeper in the atmosphere due to stronger scattering

See the report for full details and plots.

## References

Reed, L. L., J. H. Jain, W. L. Grantham, and J. D. Neff. "Radiative transfer models for geophysical applications." Applied Optics 15, no. 7 (1976): 1853-1860.

## Acknowledgements

This analysis was completed as part of individual project for ATMS 532 at University of Washington.
