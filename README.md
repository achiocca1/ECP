# ECP
```diff
- Be sure you are running the latest release available (v2.0.0)!
```
This repository hosts the implementation of the **Effective Critical Plane (ECP)** method—a novel approach for assessing fatigue strength in structural components under complex loading conditions. By averaging the stress–strain field over a small, material-dependent volume, the ECP method enhances traditional critical plane analysis. This results in improved accuracy when predicting fatigue life, especially in areas with stress concentrations such as notches and welds.

The repository includes source code, validation data, and example cases for researchers and engineers working on fatigue analysis.

The mathematical procedure is described in the article:

**A. Chiocca, F. Frendo, "Fatigue assessment of structural components through the Effective Critical Plane factor", International Journal of Fatigue, 2024, doi.org/10.1016/j.ijfatigue.2024.108565.**

---

## Repository Structure

The repository is organized into three main folders:

### _1_FEM_PostProcessing
- **ANSYS_Post_Process_Results.mac**  
  A macro for post-processing ANSYS results to prepare input data for the fatigue analysis.

### _2_CriticalPlaneEvaluation
Contains the files required to compute the traditional critical plane parameters:
- **COORD.csv**  
  Nodal coordinates file.
- **RESULTS.csv**  
  ANSYS simulation results.
- **FS_OPT.m**  
  MATLAB function to compute Fatemi–Socie parameters.
- **SWT_OPT.m**  
  MATLAB function to compute Smith–Watson–Topper (SWT) parameters.
- **MAIN_FS.m**  
  MATLAB script for processing Fatemi–Socie (FS) based calculations.
- **MAIN_SWT.m**  
  MATLAB script for processing SWT based calculations.

### _3_ECP
Contains the implementation of the ECP method and related post-processing routines:
- **MAIN_ECP.m**  
  MATLAB script that performs calibration, fatigue life prediction, and prevision data processing using the ECP approach.  
  **Note:** The script reads the polynomial fit information from a calibration file (see below) and scales the ECP factor based on the load value (third column). The idea is to run `MAIN_ECP.m` once to determine the optimal control radius and then re-run the analysis using that optimal value.
- **FatigueScatter.m**  
  MATLAB function for performing fatigue scatter calculations.

Inside the _3_ECP folder, there are also four subfolders:
- **CALIBRATION**  
  Contains calibration data files.  
  **File naming:** Files should be structured as `1_SPECIMEN_NAME_LOAD_RATIO.txt` (e.g., `_R_-1.txt`). These files contain the polynomial fit data derived from the _2_CriticalPlaneEvaluation results.
- **CORRECTOR**  
  Contains correction files, if necessary, for refining the calibration data. Be aware that the optimized method for evaluating critical plane factors, as it is currently implemented, computes amplitudes (e.g., shear strain amplitude for FS by means of FS_OPT.m) via the Minimum Circumscribed Circle (MCC) method). If I wanted to derive the amplitude via another method, e.g. Minimum Circumscribed Ellipse (MCE) then a corrective file needs to be implemented which contains a factor that increases or decreases the amplitude of the reference stress or strain with respect to the MCC.  
  **File naming:** Files should be structured as `3_Corr_SPECIMEN_NAME_LOAD_RATIO.txt`. The file must contain a single float.
- **PREVISION**  
  Contains experimental fatigue information files.  
  **File naming:** Files should be structured as `2_Exp_SPECIMEN_NAME_LOAD_RATIO.txt` (e.g., `_R_-1.txt`). The first column represents the number of cycles to fracture, the second column contains nominal normal and shear stress, and the third column provides the load or moment used during the test.
- **RESULTS**  
  This folder is used to store the output data, figures, and result tables produced during the ECP analysis.

---

## Getting Started

1. **Pre-Processing:**  
   Use the macro in the _1_FEM_PostProcessing folder (`ANSYS_Post_Process_Results.mac`) to process your ANSYS simulation data.

2. **Critical Plane Evaluation:**  
   Run the MATLAB scripts (`MAIN_FS.m` or `MAIN_SWT.m`) in the _2_CriticalPlaneEvaluation folder to compute the FS or SWT critical plane parameters.

3. **ECP Analysis:**  
   - Place your calibration and experimental data files in the corresponding subfolders inside _3_ECP following the prescribed naming conventions.
   - Run `MAIN_ECP.m` to perform the calibration and fatigue life prediction.
   - For optimal results, run `MAIN_ECP.m` twice:  
     - **First run:** Determine the optimal control radius from the calibration data.  
     - **Second run:** Use the optimal control radius (the value obtained in the first run) to scale the ECP factor based on the load value.

---

## Bug Reports & Contributions

Bug reports, suggestions, and contributions are welcome! This software is regularly maintained. If you encounter any issues or have ideas for improvements, please open an issue or submit a pull request.

---

## Contact

For further assistance, please contact Andrea Chiocca at [andrea.chiocca@unipi.it](mailto:andrea.chiocca@unipi.it).
