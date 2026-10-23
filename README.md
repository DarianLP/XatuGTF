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

In order for **XatuGTF** to run calculations for a given material, it needs two files: 
- **.outp file** — Properties output file of CRYSTAL.
- **basis.txt file** — Basis used in the properties calculation of CRYSTAL, with the same format as in the .d12 file.

The input properties file (.d3) must contain the following keywords:
```bash
BASISSET
2
60 nR
64 nR
END
```
This will print in the .outp file the Fock and overlap matrix in real space for a number nR of Bravais lattice vectors, i.e. $H(\mathbf{R}), S(\mathbf{R})$. Make sure to substitue in the .d3 file nR by the number of cells you wish to take into account. Once that is set, go into main and go to the file whose name indicates the calculation you want to perform. You will find a section with the input parameters, which you must choose. In order to compile the changes you made main, you must run from the XatuGTF directory
```bash
source xatu_config
make xatuGTF_compute_insert_specific_calculation
```
Once that is set, from the XatuGTF directory run the following
```bash
mpirun -np insert_ncores bin/xatuGTF_insert_specific_calculation
```
Once the calculation has finished, you will find the results in the file Results. 
