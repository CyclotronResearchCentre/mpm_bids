# MPM / QSM processing

This package preprocesses EPI based MPM/QSM data, performs an hMRI based MPM analysis and prepares the QSM analysis (not included, ask Monica for information)

## Installation

To install the MPMQSM_preproc use either

- *python setup.py develop*
- *python setup.py install*

## BIDS data structure

This package uses the BIDS data format. Data has to be organized before in the following manner:

```
├── MT.mat
├── PD.mat
├── cest
├── fmap
│   ├── DZNE_subj-02_ses-01_fmap-b1-con.json
│   ├── DZNE_subj-02_ses-01_fmap-b1-con.nii
│   ├── DZNE_subj-02_ses-01_fmap-b1.json
│   └── DZNE_subj-02_ses-01_fmap-b1.nii
├── magn
├── mpm
│   ├── DZNE_subj-02_ses-01_MTw_e1.json
│   ├── DZNE_subj-02_ses-01_MTw_e1.nii
│   ├── DZNE_subj-02_ses-01_MTw_e1_ph.json
│   ├── DZNE_subj-02_ses-01_MTw_e1_ph.nii
│   ├── DZNE_subj-02_ses-01_MTw_e2.json
│   ├── DZNE_subj-02_ses-01_MTw_e2.nii
│   ├── DZNE_subj-02_ses-01_MTw_e2_ph.json
│   ├── DZNE_subj-02_ses-01_MTw_e2_ph.nii
│   ├── DZNE_subj-02_ses-01_MTw_e3.json
│   ├── DZNE_subj-02_ses-01_MTw_e3.nii
│   ├── DZNE_subj-02_ses-01_MTw_e3_ph.json
│   ├── DZNE_subj-02_ses-01_MTw_e3_ph.nii
│   ├── DZNE_subj-02_ses-01_MTw_e4.json
│   ├── DZNE_subj-02_ses-01_MTw_e4.nii
│   ├── DZNE_subj-02_ses-01_MTw_e4_ph.json
│   ├── DZNE_subj-02_ses-01_MTw_e4_ph.nii
│   ├── DZNE_subj-02_ses-01_PD_e1.json
│   ├── DZNE_subj-02_ses-01_PD_e1.nii
│   ├── DZNE_subj-02_ses-01_PD_e1_ph.json
│   ├── DZNE_subj-02_ses-01_PD_e1_ph.nii
│   ├── DZNE_subj-02_ses-01_PD_e2.json
│   ├── DZNE_subj-02_ses-01_PD_e2.nii
│   ├── DZNE_subj-02_ses-01_PD_e2_ph.json
│   ├── DZNE_subj-02_ses-01_PD_e2_ph.nii
│   ├── DZNE_subj-02_ses-01_PD_e3.json
│   ├── DZNE_subj-02_ses-01_PD_e3.nii
│   ├── DZNE_subj-02_ses-01_PD_e3_ph.json
│   ├── DZNE_subj-02_ses-01_PD_e3_ph.nii
│   ├── DZNE_subj-02_ses-01_PD_e4.json
│   ├── DZNE_subj-02_ses-01_PD_e4.nii
│   ├── DZNE_subj-02_ses-01_PD_e4_ph.json
│   ├── DZNE_subj-02_ses-01_PD_e4_ph.nii
│   ├── DZNE_subj-02_ses-01_T1_e1.json
│   ├── DZNE_subj-02_ses-01_T1_e1.nii
│   ├── DZNE_subj-02_ses-01_T1_e1_ph.json
│   ├── DZNE_subj-02_ses-01_T1_e1_ph.nii
│   ├── DZNE_subj-02_ses-01_T1_e2.json
│   ├── DZNE_subj-02_ses-01_T1_e2.nii
│   ├── DZNE_subj-02_ses-01_T1_e2_ph.json
│   ├── DZNE_subj-02_ses-01_T1_e2_ph.nii
│   ├── DZNE_subj-02_ses-01_T1_e3.json
│   ├── DZNE_subj-02_ses-01_T1_e3.nii
│   ├── DZNE_subj-02_ses-01_T1_e3_ph.json
│   ├── DZNE_subj-02_ses-01_T1_e3_ph.nii
│   ├── DZNE_subj-02_ses-01_T1_e4.json
│   ├── DZNE_subj-02_ses-01_T1_e4.nii
│   ├── DZNE_subj-02_ses-01_T1_e4_ph.json
│   └── DZNE_subj-02_ses-01_T1_e4_ph.nii
└── pulses
```

this includes 
- fmap-b1 relative flipangle map of the ptx pulse (or CP, if no ptx is measured)
- fmap-b1-con maximal contrast from the $B_1$ maps (used for SPM registration)
- mpm/ $T_1$, MT and PD weighted, multi echo data


## Usage

*usage: mpmqsm_preproc [-h] --path PATH [--den DEN] --site SITE --sub SUB --ses SES [--ptx]*

```
  --path PATH  path to MPM/QSM data, folder, where the BIDS structure
               originates (default: None)
  --den DEN    boolean, shall the raw data be denoised? default = true
               (default: True)
  --site SITE  exact name of site (eg DZNE) (default: None)
  --sub SUB    exact name of subject (eg subj-01) (default: None)
  --ses SES    exact name of session (eg ses-01) (default: None)
  --ptx        create log file? default true (default: True)
```

## Output

MPM output can be found in 

*derivates/SITE/SUBJECT/SESSION/mpm/maps/Results*

QSM input (for Monica's analysis)

*derivates/SITE/SUBJECT/SESSION/qsm/ROCombine*
