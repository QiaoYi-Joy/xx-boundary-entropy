# XX-chain boundary entropy

MATLAB numerical programs accompanying the manuscript:

**Explicit Boundary-Entropy Deficits in the XX Spin Chain with Arbitrary Boundary Fields**

## Overview

This repository provides numerical programs for studying the open
spin-1/2 XX chain with general Hermitian boundary fields.

The work connects the zeros of transfer-matrix eigenvalues to
finite-size spectra and boundary thermodynamics. Its main results
include an explicit boundary free energy, a low-temperature
boundary-entropy counting law, and the temperature window in which
the corresponding contribution can be resolved in finite chains.

## Requirements

MATLAB is required. Additional toolbox requirements depend on the
functions used by the selected program.

## Usage

1. Download the repository using **Code → Download ZIP** and extract it.
2. Open the extracted folder in MATLAB, preserving the directory structure.
3. Select the calculation or plotting script relevant to the quantity
   of interest and inspect its model parameters before running it.
4. Start with a small chain length when checking finite-size calculations.

Large-chain calculations may require substantially more computation time.

## Conventions

Unless explicitly stated otherwise, the exchange coupling and
Boltzmann constant are set to J = k_B = 1.

- N: number of spin sites.
- T: temperature.
- ND: nondiagonal boundary with a nonzero transverse field.
- D: diagonal boundary with zero transverse field.
- F: free boundary with no boundary field.

Definitions of observables and analytical expressions are given in
the manuscript and its Supplemental Material.

## Citation

If you use these programs, please cite the accompanying manuscript
and identify the repository version or commit used in your calculations.

Publication details will be added when available.
