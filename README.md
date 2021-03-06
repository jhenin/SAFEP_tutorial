# SAFEP_tutorial
This repository contains all the files needed to calculate the binding free energy of phenol to lysozyme following the SAFEP Tutorial.

Structure of supplementary files:
```
Supplementary_Files
|---common
|   |---CHARMM36m_FF
|   |---structures
|       | PHEN_lyso_complex.pdb
|       | PHEN_lyso_complex.psf
|
|---Step0_unbiased_simulations
|   | PHEN_lyso_dry.pdb #or pdbid
|   | run_unbiased.namd
|
|---Step1_DefineTheBoundState
|   |---inputs
|   |   | unbiased_sample.dcd #downsampled trajectory
|   |    
|   |---sample_results
|   |   | DBC_histogram_sample.pdf
|   |   | DBC_template.colvars
|   |
|   |---outputs
|
|---Step2_alchemy_site
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
|
|---Step3_restraint_perturbation
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
|---Step4_alchemy_bulk
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
    | AFEP_parse.py
    | SAFEP_Analysis.ipynb
    | TI_Calculations.py
```
