# SAFEP_tutorial
This repository contains all the files needed to calculate the binding free energy of phenol to lysozyme following the SAFEP Tutorial.

Structure of supplementary files:
```
Supp-Files
| prot-ref.pdb #protein reference structure
| rest-ref.pdb #ligand reference structure
| solv-prot.pdb #solvated protein with ligand coordinated
| solv-prot.psf #solvated protein with ligand structure
|
|---AFEP-Bound-Decoupling
|    | DBC-Restraint.colvars #DBC colvars
|    | run.namd #run the decoupling
|    | solv-prot-charmm.fep #?
|    | 
|    |----equ
|
|---AFEP-Bulk
|---DBC-Unbiased
|---DBC_TI
|---Force-field-parameters
|---Python-scripts
|---Jupyter-notebook
|---common
|
```

Proposed structure:
```
Supplementary_Files
|---common
|   |---CHARMM36m_FF
|   |---Structures
|       | prot_solv.pdb
|       | prot_solv.psf
|
|---Step1_DefineTheBoundState
|   |---inputs
|   |   | run_unbiased.namd
|   |    
|   |---sample_outputs
|   |   | unbiased_sample.dcd #downsampled trajectory
|   |   | DBC_sample.colvars
|   |   | DBC_histogram_sample.pdf
|   |
|   |---outputs
|
|---Step2_FEP_dG_site
|   |---inputs
|   |   | DBC_sample.colvars
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
|---Step3_TI_dG_DBC
|   |---inputs
|   |   | DBC_sample.colvars
|   |   | run_siteFEP.namd
|   |   | prot_solv.fep.pdb
|   |   | prot_solv.ref.pdb
|   |    
|   |---sample_outputs
|   |   | siteFEP.fepout
|   |
|   |---outputs
|
|---Step4_FEP_dG_bulk
|   |---inputs
|   |   | DBC_sample.colvars
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
