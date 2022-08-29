# SAFEP_tutorial
This repository contains all the files needed to calculate the binding free energy of phenol to lysozyme following the SAFEP Tutorial.

Structure of supplementary files:
```

|---common
|   |---CHARMM36m_FF
|   |---Python_Modules
|   |---structures
|       | phenol_lysozyme.psf
|       | phenol_lysozyme.pdb
|       | phenol_water.psf
|       | phenol_water.pdb
|       | phenol_gas_phase.psf
|       | phenol_gas_phase.pdb
|   | fep.tcl
|
|---StepA_create_DBC
|   |---inputs
|       | unbiased_01.dcd #downsampled trajectory
|   |---sample_outputs
|       | DBC_restraint.colvars
|   |---Step0_unbiased_simulations
|   |   |---inputs
|           | run_unbiased.namd
|
|---StepB_alchemy_site
|   |---inputs
|   |   | DBC_template.colvars
|   |   | run_siteFEP.namd
|   |    
|   |---sample_outputs
|   |   | siteFEP.fepout
|   |
|   |---outputs
|
|
|---StepC_restraint_perturbation
|   |---inputs
|   |   | DBC_template.colvars
|   |   | restraint_perturbation.namd
|   |   | gas_phase_PHEN.fep.pdb
|   |    
|   |---sample_outputs
|   |   | restraint_perturbation.colvars.traj
|   |
|   |---outputs
|
|---StepD_alchemy_bulk
|   |---inputs
|   |   | DBC_template.colvars
|   |   | run_siteFEP.namd
|   |   | prot_solv.fep.pdb
|   |   | prot_solv.ref.pdb
|   |    
|   |---sample_outputs
|   |   | siteFEP.fepout
|   |
|   |---outputs
|
|---Jupyter_Notebook
|    | SAFEP_Analysis.py
|    | TI_Calculations.py
|    
|---Text_Source
|   |
|   | .tex files
|   | figures
|
| Tutorial.pdf
| SAFEP 
```

