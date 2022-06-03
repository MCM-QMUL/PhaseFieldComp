# PhaseFieldComp
**Open-source Code for Phase Field Damage Modelling of Composite Materials**

![alt text](https://github.com/MCM-QMUL/PhaseFieldComp/blob/32bcf04f88e2abe90b0683cafd7f88799b3763bb/docs/PFM_microscale.gif)

PhaseFieldComp is an open-source framework capable of carrying out phase field (PF) fracture simulations. The formulation combines two fracture models. The phase field fracture method, capable of capturing arbitrary crack trajectories, is used to model crack initiation and growth along the matrix and the fibres. Furthermore, fibre-matrix debonding is simulated using a cohesive zone model. The framework is implemented in the commercial finite element package ABAQUS by user subroutines.

## Features
- Beginner friendly framework.
- Can be adapted to model other microscale fracture in composite materials.
- ABAQUS, user-defined Fortran subrotuines, Python.
- Suitable for parrallel computing.
- Suitable for Windows, Linux, MacOX systems.

## Limitations
- Linear elastic fracture
- Models that consider plasticicity will be published in the future

## Requirements or Installation 

Follow this [link](https://bibekanandadatta.com/link-intel-and-vs-abaqus-2020/) to install and link the following three software.
1. Install ABAQUS
2. Visual Studio
3. Fortran complier 

4. Then, create your CAE model and generate an original input file. 
5. Prepare your ABAQUS input (.inp) files using the Python scripts in the "./script" folder. The input file examples (original and ready-for-PFM) are in the script folder. Run the python script will generate an input file suitable for PF model. 
6. Submit your job via ABAQUS command line: abaqus analysis job=YOURFILE user=subroutine.for standard_parallel=solver double=both cpus=x (x is the number of CPU cores).
7. Then you can visuaslise the results via ABAQUS. Check the documentation and paper for more informaiton. 

## Reference
If using this code for research or industrial purposes, please cite:
[1] Tan, W. and Martínez-Pañeda, E., 2021. Phase field predictions of microscopic fracture and R-curve behaviour of fibre-reinforced composites. Composites Science and Technology, 202, p.108539. doi: https://doi.org/10.1016/j.compscitech.2020.108539
[2] Tan, W. and Martínez-Pañeda, E., 2022. Phase field fracture predictions of microscopic bridging behaviour of composite materials. Composite Structures, 286, p.115242. doi: https://doi.org/10.1016/j.compstruct.2022.115242

## License
MIT
