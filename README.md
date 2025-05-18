# PaCo in SageMath

This repository provides a proof-of-concept implementation in SageMath
of the PaCo bootstrapping procedure [1].

## Table of Contents

- [Installation](#installation)
- [Repository Structure](#repository-structure)
- [Usage](#usage)
- [References](#references)

## Installation

1. **Prerequisites**:
Ensure you have SageMath 10.0 or higher installed.
You can download it from the
[SageMath website](https://www.sagemath.org/download.html).
   
2. **Clone the repository**:
   ```bash
   git clone https://github.com/se-tim/PaCo-Implementation.git
   cd PaCo-Implementation
   ```

3. **Initialize the submodules**:
Ensure all submodules are correctly initialized,
as explained here:
[Git submodule inside of a submodule (nested submodules)](https://stackoverflow.com/q/1535524).
Therefore, run the following command inside the cloned repository:
   ```bash
   git submodule update --init --recursive
   ```

## Repository Structure

- **`paco_package/`**:
The main package containing all implementations
related to CKKS and the PaCo procedure:
   - **`ckks_in_sagemath/`**:
   A submodule.
   It implements the main CKKS functionalities,
   including the original bootstrapping procedure.
   Make sure that this is initialized correctly,
   as described in the section [Installation](#installation).
  
  - **`ckks_x.py`**:
   Implements the functionalities required for the PaCo procedure,
   building on the core CKKS operations.

  - **`ext_bit_rev.py`**:
  Implements functions related to extended bit-reversing.

  - **`fast_dft_with_ext_bit_rev.py`**:
  Provides the matrices required for the partial CoeffToSlot
  and the conventional SlotToCoeff operations in the PaCo procedure.

- **`benchmarks.py`**
A script evaluating the performance of the PaCo procedure
against the original bootstrapping procedure.

## Usage

The `benchmarks.py` script compares the performance of the PaCo procedure
with the original CKKS bootstrapping procedure:

1. **Parameter selection**:  
   Choose between predefined parameter sets (from the paper [1])
   or input custom parameters.
   Further,
   specify how many test runs to perform for each configuration,
   and select either sequential or parallel bootstrapping.
   If parallel bootstrapping is chosen,
   you will also be prompted to enter the number of simulated processors
   to use for the PaCo procedure.

2. **Performance output**:
   A summary table comparing the timings
   of the two bootstrapping methods is output.
   For the sequential bootstrapping, 
   the average precision is also computed.

To run the script,
make sure you are in the root directory of the repository,
then execute:
```bash
sage benchmarks.py
```
Results will be printed directly to the terminal.

## References

1. (Add later)