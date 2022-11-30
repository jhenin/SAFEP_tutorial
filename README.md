# SAFEP_tutorial
This repository contains all the files needed to calculate the binding free energy of phenol to lysozyme following the SAFEP Tutorial.

Structure of supplementary files:
```
├── Binding_Tutorial.pdf
├── common
│   ├── CHARMM36m_FF
│   │   ├── par_all36_carb.prm
│   │   ├── par_all36_cgenff.prm
│   │   ├── par_all36_lipid.prm
│   │   ├── par_all36_na.prm
│   │   ├── par_all36_prot.prm
│   │   ├── toppar_all36_lipid_cationpi_wyf.str
│   │   └── toppar_water_ions_namd.str
│   ├── common_config.namd
│   ├── fep.tcl
│   ├── protein_tilt.colvars
│   ├── structures
│   │   ├── phenol_gas_phase.pdb
│   │   ├── phenol_gas_phase.psf
│   │   ├── phenol_lysozyme.pdb
│   │   ├── phenol_lysozyme.psf
│   │   ├── phenol_water.pdb
│   │   └── phenol_water.psf
│   ├── TI.tcl
│   └── TI_template.colvars
├── file_structure.txt
├── README.md
├── SAFEP_Tutorial_Notebook.ipynb
├── stepA_create_DBC
│   ├── inputs
│   │   ├── unbiased_sample.dcd
│   │   └── vmdscene.tga
│   ├── sample_outputs
│   │   └── DBC_restraint.colvars
│   └── step0_unbiased_simulation
│       ├── inputs
│       │   └── run_unbiased.namd
│       └── outputs
├── stepB_alchemy_site
│   ├── inputs
│   │   └── run.namd
│   ├── outputs
│   └── sample_outputs
│       ├── alchemy_site.fepout
│       ├── alchemy_site.pdb
│       ├── bound_convergence.pdf
│       └── bound_generalFigures.pdf
├── stepC_restraint_perturbation
│   ├── inputs
│   │   └── run.namd
│   ├── outputs
│   └── sample_outputs
│       ├── DBC_restraint_RFEP.colvars
│       ├── RFEP.colvars.traj
│       └── TI_general.pdf
├── stepD_alchemy_bulk
│   ├── inputs
│   │   └── run.namd
│   ├── outputs
│   └── sample_outputs
│       ├── alchemy_bulk.fepout
│       ├── alchemy_bulk.pdb
│       ├── bulk_convergence.pdf
│       ├── bulk_convergence.svg
│       └── bulk_generalFigures.pdf
├── text_src
│   ├── Figures
│   │   ├── AFEP2.png
│   │   ├── AFEP2-prob.png
│   │   ├── AFEP2-solution.png
│   │   ├── AFEP-decoupling-lambda.png
│   │   ├── AFEP-decoupling-sum.png
│   │   ├── DBC.png
│   │   ├── DBCsym.png
│   │   ├── example_symmetry_labels.png
│   │   ├── histogram_nosym.png
│   │   ├── histogram.png
│   │   ├── HSELEU.jpg
│   │   ├── ipynb-cell4.png
│   │   ├── ipynb.png
│   │   ├── lyso5-new.jpg
│   │   ├── phenol-permutation.jpg
│   │   ├── phenol-permutation-new.jpg
│   │   ├── probability.png
│   │   ├── RFEP-log.png
│   │   ├── RFEP.png
│   │   ├── scheme.jpg
│   │   ├── scheme.vsdx
│   │   ├── thermo-cycle.jpg
│   │   └── Thermo-cycle.vsdx
│   ├── livecoms.cls
│   ├── livecoms-template-tutorials.tex
│   ├── Ref.bib
│   ├── tutorial.tex
│   └── vancouver-livecoms.bst
└── titration_curve.pdf

23 directories, 70 files
```

