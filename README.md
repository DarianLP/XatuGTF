# XatuGTF

**XatuGTF** (*eXcitons from ATomistic calcUlations with a Gaussian Type Function basis*) is a scientific software package designed to solve the **Bethe–Salpeter equation (BSE)** for a wide range of materials.

Starting from electronic band structures obtained from density functional theory (DFT) calculations using a Gaussian-type orbital basis, XatuGTF constructs and diagonalizes the Bethe–Salpeter Hamiltonian to obtain the exciton spectrum. The resulting excitonic states can then be post-processed to characterize their properties and investigate the optical response of the material.

## Installation

### 1. Install the dependencies

XatuGTF requires the following numerical libraries:

- **OpenBLAS** — optimized linear algebra routines.
- **LAPACK** — numerical linear algebra routines.
- **ARPACK** — eigenvalue and eigenvector computation.
- **Armadillo** — a C++ linear algebra library.
- **Open MPI** — an implementation of the Message Passing Interface (MPI) for parallel computing.

On Ubuntu, install the dependencies by running:

```bash
sudo apt update
sudo apt install libopenblas-dev liblapack-dev libarpack2-dev libarmadillo-dev
sudo apt install openmpi-bin libopenmpi-dev
```

### 2. Compile the code

Navigate to the XatuGTF root directory:

```bash
cd /path/to/XatuGTF
```

Build by running:

```bash
make build
make xatu
```

Once compilation completes successfully, XatuGTF is ready to use.

## Usage

Instructions for running calculations and post-processing the results will be provided here.