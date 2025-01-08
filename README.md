# Simulation of radar-wave propagation in disaster sites for drone-based life detection
 GitHub repository for the source files used for the creation of Dominik Spale's master thesis with the title "Simulation of radar-wave propagation in disaster sites for drone-based life detection". Includes all necessary files for gprMax simulations as well as MATLAB data processing.


## Installation
The simulations were performed in [gprMax](https://www.gprmax.com/), an open-source FDTD simulation program. To install gprMax, please follow the instructions given at their respective [GitHub page](https://github.com/gprmax/gprMax).

To install MATLAB, please follow the instructions on the [MathWorks website](https://de.mathworks.com/products/matlab/getting-started.html).
## Usage
To perform a simulation in gprMax, execute the respective scripts with gprMax:
- For monochromatic continuous wave radar simulations (MCCW), refer to [Monochromatic Continuous Wave](<./Monochromatic Continuous Wave>).
- For FCCW bioradar simulations with the superposition excitation approach, refer to [FCCW - Superposition excitation](<./FCCW - Superposition excitation>).
- For FCCW bioradar simulations with the single excitation approach, refer to [FCCW - Single excitation](<./FCCW - Single excitation>).

If the simulation is run in a HPC environment (as recommended), examples for shell scripts for the bwHPC cluster are given in [Shell scripts](<./Shell scripts>). Adjust the respective parameters according to your folder structure.

The MATLAB processing files are given in [MATLAB files](<./MATLAB files>). Additionally, MATLAB scripts for amplitude and phase determination via linear regression and Levenberg-Marquardt optimization for circle fitting are included. To perform FCCW simulations with the single excitation approach, the excitation signal needs to generated in a .txt-file first via [FCCW_generation.m](<./MATLAB scripts/FCCW_generation.m>).

To process data, place the respective evaluation script into the same folder as the output files of the gprMax simulations. Adjust image_path (where images are stored) and plotting parameters according to your needs.

## Feedback

For feedback, bug reports, etc., please reach out to me via dominik@spale.imtek-uni-freiburg.de. Thanks!


## 🔗 Links
[![linkedin](https://img.shields.io/badge/linkedin-0A66C2?style=for-the-badge&logo=linkedin&logoColor=white)](www.linkedin.com/in/dominik-martin-spale-51bb0a202/)