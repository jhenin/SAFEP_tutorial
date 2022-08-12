<<<<<<< HEAD
# SAFEP_tutorial
This repository contains all the files needed to calculate the binding free energy of phenol to lysozyme following the SAFEP Tutorial.

Structure of supplementary files:
```

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
|---Step1_create_DBC
|   |---inputs
|   |   | unbiased_sample.dcd #downsampled trajectory
|   |    
|   |---sample_outputs
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
=======
# SAFEP_tutorial
This repository contains all the files needed to calculate the binding free energy of phenol to lysozyme following the SAFEP Tutorial.

Structure of supplementary files:
```
├── Binding_Tutorial.pdf
├── Figures
│   ├── AFEP2.png
│   ├── AFEP2-prob.png
│   ├── AFEP2-solution.png
│   ├── AFEP-decoupling-lambda.png
│   ├── AFEP-decoupling-sum.png
│   ├── DBC.png
│   ├── DBCsym.png
│   ├── histogram_nosym.png
│   ├── histogram.png
│   ├── HSELEU.jpg
│   ├── ipynb-cell4.png
│   ├── ipynb.png
│   ├── lyso5-new.jpg
│   ├── phenol-permutation.jpg
│   ├── phenol-permutation-new.jpg
│   ├── probability.png
│   ├── RFEP-log.png
│   ├── RFEP.png
│   ├── scheme.jpg
│   ├── scheme.vsdx
│   ├── thermo-cycle.jpg
│   └── Thermo-cycle.vsdx
├── livecoms.cls
├── livecoms-template-tutorials.tex
├── README.md
├── Ref.bib
├── Supp-Files
│   ├── common
│   │   ├── CHARMM36m_FF
│   │   │   ├── par_all36_carb.prm
│   │   │   ├── par_all36_cgenff.prm
│   │   │   ├── par_all36_lipid.prm
│   │   │   ├── par_all36_na.prm
│   │   │   ├── par_all36_prot.prm
│   │   │   └── toppar_water_ions_namd.str
│   │   ├── fep.tcl
│   │   ├── Python_Modules
│   │   │   ├── AFEP_parse.py
│   │   │   ├── SAFEP_Analysis.py
│   │   │   └── TI_Calculations.py
│   │   └── structures
│   │       ├── PHEN_lyso_complex.pdb
│   │       ├── PHEN_lyso_complex.psf
│   │       ├── PHEN.pdb
│   │       └── PHEN.psf
│   ├── SAFEP_Tutorial_Notebook.ipynb
│   ├── Step0_unbiased_simulation
│   │   └── inputs
│   │       └── run_unbiased.namd
│   ├── Step1_create_DBC
│   │   ├── inputs
│   │   │   └── unbiased-01.dcd
│   │   └── sample_outputs
│   │       └── DBC_restraint.colvars
│   ├── Step2_alchemy_site
│   │   ├── inputs
│   │   │   ├── fep.tcl
│   │   │   └── run.namd
│   │   ├── output
│   │   └── sample_outputs
│   │       ├── alchemy_site.fep
│   │       ├── bound_convergence.pdf
│   │       └── bound_generalFigures.pdf
│   ├── Step3_restraint_perturbation
│   │   ├── DBC-Restraint-RFEP.colvars
│   │   ├── DBC-Restraint-RFEP_new.colvars
│   │   ├── ligand-ref.pdb
│   │   ├── output
│   │   ├── PHEN.pdb
│   │   ├── PHEN.psf
│   │   ├── run.namd
│   │   └── sample_output
│   │       ├── RFEP.colvars.traj
│   │       ├── RFEP.dat
│   │       ├── RFEP.dcd
│   │       ├── RFEP.log
│   │       ├── TI_general.pdf
│   │       └── TI_general.svg
│   ├── Step4_alchemy_bulk
│   │   ├── equ_output
│   │   │   ├── AFEP-Hyd-01.restart.coor
│   │   │   ├── AFEP-Hyd-01.restart.vel
│   │   │   ├── AFEP-Hyd-01.restart.xsc
│   │   │   ├── AFEP-Hyd.log
│   │   │   └── equ.namd
│   │   ├── output
│   │   ├── PHEN-solv.fep
│   │   ├── PHEN-solv.pdb
│   │   ├── PHEN-solv.psf
│   │   ├── run.namd
│   │   └── sample_output
│   │       ├── AFEP-Hyd-02.coor
│   │       ├── AFEP-Hyd-02.dcd
│   │       ├── AFEP-Hyd-02.fepout
│   │       ├── AFEP-Hyd-02.restart.coor
│   │       ├── AFEP-Hyd-02.restart.vel
│   │       ├── AFEP-Hyd-02.restart.xsc
│   │       ├── AFEP-Hyd-02.vel
│   │       ├── AFEP-Hyd-02.xsc
│   │       ├── AFEP-Hyd-02.xst
│   │       ├── bulk_convergence.pdf
│   │       ├── bulk_convergence.svg
│   │       ├── bulk_generalFigures.pdf
│   │       └── VMD-ParseFEP
│   │           ├── deinterleave_idws.py
│   │           ├── free-energy.1.png
│   │           ├── free-energy.2.png
│   │           ├── ParseFEP.log
│   │           ├── probability.1.png
│   │           ├── probability.2.png
│   │           └── summary.png
│   └── titration_curve.pdf
├── tutorial.tex
└── vancouver-livecoms.bst

```
>>>>>>> dd06a8b8b62f4661bfe7aef6c140b678f4fe9641
