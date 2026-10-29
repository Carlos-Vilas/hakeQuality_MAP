# Mathematical modeling of quality degradation in hake fillets under MAP

This repository contains the MATLAB code and experimental dataset used in the study:

> **Mathematical modeling of the degradation of biochemical and subjective quality indicators in Hake fillets stored under modified atmosphere packaging**

The code implements a hierarchical framework that links storage temperature to objective fish-quality indicators through knowledge-driven dynamic models and subsequently relates those indicators to the sensory Quality Index Method (QIM) through data-driven regression models. Modified Atmosphere packaging is considered. The objective indicators considered are *Shewanella* spp., *Photobacterium phosphoreum*, total volatile basic nitrogen (TVB-N), and nucleotide degradation products. 

## Repository contents

```text
.
├── dataDrivenModels/
│   ├── kvalue/
│   │   ├── mlp_cv/
│   │   │   ├── model/
│   │   │   │   ├── defineMLP.m
│   │   │   │   └── evaluateModelPerform.m
│   │   │   ├── preprocessing/
│   │   │   │   ├── buildData.m
│   │   │   │   └── generateSyntheticData.m
│   │   │   ├── results/
│   │   │   ├── mainSelectNN.m
│   │   │   └── mainTrainNN.m
│   │   ├── plsr_lin/
│   │   ├── plsr_nl/
│   │   └── randomForest/
│   └── nucleotide/
│       ├── mlp_cv/
│       ├── plsr_lin/
│       ├── plsr_nl/
│       └── randomForest/
├── experimentalData/
│   ├── hakeQuality_expData.db
│   └── readExperimentalData.m
├── knowledgeDrivenModels/
│   ├── bads/
│   │   ├── model01/
│   │   │   ├── computeNLL_reps.m
│   │   │   ├── costWrapper_01.m
│   │   │   ├── mainPE_model01.m
│   │   │   ├── mainSim_model01.m
│   │   │   ├── mex_model01.mexa64
│   │   │   ├── PEresults_model01.mat
│   │   │   └── simulationModel01.m
│   │   ├── model02/
│   │   ├── model03/
│   │   ├── model04/
│   │   ├── model05/
│   │   ├── model06/
│   │   └── model07/
│   └── vbmc/
│       └── model07/
├── mexFiles/
└── strike-goldd-files/
```

Folders representing analogous model variants or data-driven methods follow structures similar to the expanded examples above.

## Folder details

### `experimentalData/`

- `hakeQuality_expData.db`: database containing the experimental data used by the code.
- `readExperimentalData.m`: MATLAB routine for reading and organizing the experimental data.

The experiments include storage-temperature profiles and measurements used in the modeling workflow for bacterial populations, TVB-N, nucleotide compounds, and QIM. 

### `knowledgeDrivenModels/`

Contains the dynamic models for the objective quality indicators.

#### `bads/`

Contains seven candidate model structures (`model01`-`model07`) and the files used for parameter estimation and simulation with Bayesian Adaptive Direct Search (BADS). All folders follow the same general organization:

- `mainPE_model0X.m`: main script for parameter estimation.
- `mainSim_model0X.m`: main script for model simulation.
- `simulationModel0X.m`: implementation of the dynamic model used in simulations.
- `computeNLL_reps.m`: evaluates the negative log-likelihood using replicated measurements.
- `costWrapper_0X.m`: objective-function wrapper used during optimization.
- `mex_model0X.mexa64`: compiled model binary for 64-bit Linux.
- `mex_model0X.mexw64`: compiled model binary for 64-bit Windows.
- `PEresults_model0X.mat`: stored parameter-estimation output.

#### `vbmc/`

Contains the Variational Bayesian Monte Carlo workflow used to quantify parameter uncertainty for `model07`.  Files  `simulationModel07.m` ,  `computeNLL_reps.m`, `mex_model0X.mexa64`, `mex_model0X.mexw64`, and `PEresults_model0X.mat` are the same as in the BADS folder. The new files in VBMC are:

- `mainCI_model07.m`: main script for parameter uncertainty estimation.
- `costWrapper_07_logJoint.m`: objective-function (log-posterior) wrapper used during optimization.
- `PEresults_model07_CI.mat`: stored uncertainty estimation output.

The mechanistic framework includes temperature-dependent bacterial-growth models, alternative TVB-N formulations, and nucleotide-degradation models with candidate reaction pathways and bacterial contributions. 

### `dataDrivenModels/`

Contains the regression approaches used to map objective quality indicators to QIM:

- `mlp_cv/`: multilayer perceptron with cross-validation-based model selection.
- `plsr_lin/`: linear partial least-squares regression.
- `plsr_nl/`: nonlinear partial least-squares regression using an RBF kernel.
- `randomForest/`: random-forest regression.

Two predictor formulations are evaluated:

- `kvalue/`: uses the K-index representation together with bacterial and TVB-N predictors.
- `nucleotide/`: uses the nucleotide compounds separately together with bacterial and TVB-N predictors.

In the expanded `mlp_cv/` example:

- `mainSelectNN.m`: main script for neural-network model selection. 
- `mainTrainNN.m`: main script for neural-network training.
- `model/defineMLP.m`: defines the multilayer-perceptron architecture.
- `model/evaluateModelPerform.m`: evaluates model performance using the MAE and R2.
- `preprocessing/buildData.m`: constructs the datasets (modeling and test) used by the model.
- `preprocessing/generateSyntheticData.m`: generates synthetic observations for model development.
- `results/`: destination for generated outputs.

The other folders follow structures similar to the `mlp_cv/` example.

### `mexFiles/`

Contains the compiled versions of the MEX files (Linux and Windows 64 bits) used to integrate the ordinary differential equations describing the evolution of the objective quality indicators. 

It also contains the MATLAB codes used to generate the C version of the different knowledge-driven models and to compile them (together with the CVODES codes) to generate the MEX files. To use these codes, follow the next steps

1.- Download the [AMIGO2 Toolbox](https://sites.google.com/site/amigo2toolbox/download)

2.- Unzip the `AMIGO2_R2025.zip` file

3.- Copy the content of folder `AMIGO2_R2025/Kernel/IVP_solvers/cvodes/C_src4Amigo`  into the `mexFiles/C_src4Amigo` folder


### `strike-goldd-files/`

Contains files required to evaluate structural identifiability with STRIKE-GOLDD. The steps to use these files are:

1.- Download the STRIKE-GOLDD toolbox.

2.- Copy the files `create_matFile_fishQuality.m` and  `mainCreateFQ_matFiles.m` into the `models` folder of STRIKE-GOLDD

3.- Copy the file `options_fishQuality.m` into the root path of STRIKE-GOLDD

4.- Generate the mat file used for the analysis (it contains the model/parameter/states/inputs... information). To that purpose:

- Open the script that creates the MATLAB data file (*.mat) that Strike-Goldd will use to analyze the model structural identifiability. 
	- The name of the file is *mainCreateFQ_matFiles.m*
	- Choose the model name ('M01', 'M02', etc)
	- Choose the parameter set: 'Complete' or 'Partial'.
- Run the matlab script
- Repeat these steps for the other models
  
5.- Set the options file for identifiability analysis

- open file options_fishQuality.m
- Choose the modelNumber and the paramSet acoording to your needs

6.- Run the identifiability analysis. In matlab 
	
	>> install	
   	>> STRIKE_GOLDD('options_fishQualtiy.m')

Results are saved in the results folder

## Requirements

- MATLAB R2025b.
- Required MATLAB toolboxes: Parallelization toolbox, Statistics and Machine Learning Toolbox, Deep Learning Toolbox
- BADS:  [Bayesian Adaptive Direct Search - v1.1.2](https://github.com/acerbilab/bads.git)
- VBMC: [Variational Bayesian Monte Carlo - v1.0.12](https://github.com/acerbilab/vbmc).
- STRIKE-GOLDD: [STRIKE-GOLDD - v4.3.0](https://github.com/afvillaverde/strike-goldd).
- Supported operating system: Linux/Windows.

The code uses MATLAB functionality including `plsregress`, `TreeBagger`, and neural-network training routines. Availability of these functions depends on the installed MATLAB products. 

> [!NOTE]
> Some scripts may expect a particular working directory or relative folder structure. Preserve the repository hierarchy unless the path handling in the scripts is updated accordingly.

## Experimental data

The repository includes the experimental database required by the supplied analyses. The underlying study comprises four independent refrigerated-storage experiments with hake fillets under modified atmosphere packaging. 

## Authors 

- Carlos Vilas,  Silvia Muñoz, Miguel Novoa, Juan R. Herrera, Marta Bernárdez, María M. Costa, Míriam R. García

Affiliations and contribution roles are provided in the accompanying manuscript. 

## Funding

This work acknowledges support from the Pescando 4.0, DSInVar, Programa de Ciencias Marinas (ThinkInAzul), and Xunta de Galicia IN607B 2026/03 projects. Complete funding statements and grant details are provided in the manuscript. 

## License

This project is licensed under the GNU General Public License v3.0.

This software is provided without warranty of any kind. The authors shall not be held liable for any damages arising from its use. See the **LICENSE** file for the complete terms.


