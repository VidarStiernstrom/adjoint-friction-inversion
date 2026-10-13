# Testing and verification of numerical methods in Rupture sensitivity to dynamic source parameters revealed by variational fracture mechanics and adjoint rupture dynamics
This repository contains the adjoint-based code for earthquake rupture with slip-weakening friction. The content is presented in the paper
'Rupture sensitivity to dynamic source parameters revealed by variational fracture mechanics and adjoint rupture dynamics' by Fukushima and Dunham.

Note that this code is an extension of Stiernström and Almquist (2024), who developed the adjoint-based dynamic rupture inversion code with rate-and-state friction.

Stiernström, V. & Almquist, M. Testing and Verification of Numerical Methods in Adjoint-Based Inversion for Stress and Frictional Parameters in Earthquake Modeling. (Zenodo, 2024). doi:10.5281/ZENODO.13776965.

## Getting started
To run the scripts, first add subdirectories to the Matlab path by running `addSubpaths` in the Matlab command window.

## Reproducing figures from the paper
- Inverse crime for initial stress: Run the script `runSWinversion.m`. Inverts for the initial stress tau0, using synthetic data recorded from a forward solve on the same computational grid. The figure is reproduced by running `runPlotSWinversion.m`.
- Sensitivity test for rupture time: Run the script `runSWvariation.m` and make the figure by running `runPlotSWvariation.m`. Computes the gradient of the rupture time 'Tdc', defined as the time when the slip reaches the characteristic slip distance 'Dc', from a reference forward simulation.

## Parameters / functions for slip-weakening friction
   - `optimize_antiplaneshear_sw.m` - Function used to invert for sw friction parameters using synthesized data, modified from `optimize_antiplaneshear.m`
   - `+pars/swFriction2DHorizontalFaultInversion.m` - Parameter set used for sw inversions, modified from `+pars/rsFriction2DHorizontalFaultInversion.m`

## Further details on the implementation:
   - `AntiplaneShear2DSWFrictionOpt.m` - Optimization class
   - `+elastic/AntiplaneShear2DSWFrictionFwdDiscr.m` - Discretization class of forward problem
   - `+elastic/AntiplaneShear2DSWFrictionAdjDiscr.m` - Discretization class of adjoint problem
   - `+elastic/+friction/+swinversion/` - Contains the slip-weakening friction functions (strength evolution), together with their derivatives w.r.t. the parameters.

## Other new files from the previous version
   - `+adjopt/optimize_func.m` - Modified from `+adjopt/optimize_lbfgs.m` to be generalized to other optimization methods, such as L-BFGS, BFGS, and gradient descent.
   - `+elastic/+discrs/+antiplaneshear/discrOps2Dsw.m` - Modified from `discrOps2D.m`: the variable name is changed from Psi (state variable in RSF) to s (slip).


# Original readme　: Testing and verification of numerical methods in Adjoint-based inversion for stress and frictional parameters in earthquake modeling
This repository contains Matlab scripts and functions for testing and verification of the numerical methods presented in the paper 'Adjoint-based inversion for stress and frictional parameters in earthquake modeling' by Stiernström et al.



**Note**:
The code in this repository is based on a bare copy of the Mercurial repository [elastic-optimization](https://sourceforge.net/p/elastic-optimization/code/ci/zenodo-archive/tree/) available on SourceForge. The code, subrepositories included, is reproduced here with the permission of the original authors. Please note that this repository does not contain the original Mercurial version history. Links to the Mercurial repositories and the changeset IDs are found below. For the original version controlled state, including history, we refer to the SourceForge repository.

## Before cloning
First install and initialize [Git Large File Storage](https://git-lfs.com). Large .mat-files, which are stored remotely, will then be downloaded when you clone the repository. Without git LFS installed, only a pointer is downloaded.

## Getting started
To run the scripts, first add subdirectories to the Matlab path by running `addSubpaths` in the Matlab command window.

## Reproducing figures from the paper
- Verification of discrete gradient: Run the script `gradientConvergenceFriction.m`. Performs a convergence study of the discrete gradient with displacement and velocity misfits. This is done through successive refinements of a finite difference approximation of the gradient.
- Inverse crime for direct effect: Run the script `direct_effect_inverse_crime.m`. Inverts for the direct effect parameter 'a', using synthetic data recorded from a forward solve on the same computational grid.

## Other useful scripts for testing/verification
   - `gradientTestFriction2D.m` - Compares the adjoint gradient to a single finite difference gradient evaluation
   - `optimize_antiplaneshear.m` - Function used to invert for friction parameters using synthesized data.
   - `fractal_fault_2d.m` - Runs a forward simulation plotting e.g. wavefields.
   - `fractal_fault_generate_data.m` - Script used to generate and save synthetic data from a forward solve.

## Parameters used for running tests are found in the subfolder `+pars`.
   - `+pars/rsFriction2DFractalFaultVerification.m` - Parameter set used for gradient verification
   - `+pars/rsFriction2DFractalFaultInversion.m` - Parameter set used for inversions.
   - `+pars/rsFriction2DFractalFaultForward.m` - Parameter set for high-resolution simulation. Used to generate high-resolution synthetic data.

## Further details on the implementation:
   - `AntiplaneShear2DRSFrictionOpt.m` - Optimization class
   - `+elastic/AntiplaneShear2DRSFrictionFwdDiscr.m` - Discretization class of forward problem
   - `+elastic/AntiplaneShear2DRSFrictionAdjDiscr.m` - Discretization class of adjoint problem
   - `+elastic/+friction/+inversion/` - Contains the friction functions (friction law, state evolution). These together with derivative of the functions w.r.t. their parameters are generated through the function `generateFunctions.m` in the same directory.

## Original Mercurial repositories on SourceForge
Links to the original Mercurial repositories and the changeset IDs used to create this repository.
   - [elastic-optimization](https://sourceforge.net/projects/elastic-optimization/): `7cc3931c602f`
   - [adjoint-optimization-package](https://sourceforge.net/projects/adjoint-optimization-package/): `5b69109e0a44`
   - [elastic-package](https://sourceforge.net/projects/elastic-package/): `10cd35e8227d`
   - [sbplib](https://sourceforge.net/projects/sbplib/): `5bffccf78673`
