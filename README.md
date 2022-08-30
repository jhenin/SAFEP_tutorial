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
|       | unbiased_sample.dcd #downsampled trajectory
|   |---sample_outputs
|       | DBC_restraint.colvars
|   |---Step0_unbiased_simulations
|   |   |---inputs
|           | run_unbiased.namd
|
|---StepB_alchemy_site
|   |---inputs
|   |   | run.namd
|   |    
|   |---sample_outputs
|   |   | alchemy_site.pdb
|   |   | alchemy_site.fepout
|   |   | bound_convergence.pdf
|   |   | bound_generalFigures.pdf
|   |
|   |---outputs
|
|
|---StepC_restraint_perturbation
|   |---inputs
|   |   | run.namd
|   |    
|   |---sample_outputs
|   |   | DBC_template.colvars
|   |   | restraint_perturbation.colvars.traj
|   |
|   |---outputs
|
|---StepD_alchemy_bulk
|   |---inputs
|   |   | run.namd
|   |    
|   |---sample_outputs
|   |   | alchemy_bulk.pdb
|   |   | alchemy_bulk.fepout
|   |   | bulk_convergence.pdf
|   |   | bulk_generalFigures.pdf
|   |
|   |---outputs
|
|---Text_src
|   |
|   | .tex files
|   | figures
|
| Binding_Tutorial.pdf
| SAFEP_Tutorial_Notebook.ipynb
| titration_curve.pdf
| README.md
```

