import os
import numpy as np
import nibabel as nib

from scipy.ndimage import gaussian_filter

"""
This scripts creates readout combined images from the raw MPM/QSM acquisitions for ALL sites/subjects/sessions/acquisition
It saves these data in the format that is needed for MPM and QSM in derivatives/subject/site/session/mpm and derivatives/subject/site/session/mpm, respectively

change path to BIDS in the very last line of this file
"""

def combine_mpm(in_folder, out_folder):
    """
    1. find all pairs of magnitide/phase images in input folder
    2. call combine_complex for all pairs
        a. combine readout directions
        b. save readout combined images to output folder
    """
    files = [f for f in os.listdir(in_folder) if ((".nii" in f) and not ("ph.nii" in f))]
    for f in files:
        combine_complex(f, in_folder, out_folder)
        print("RO combined %s" %os.path.join(in_folder,f))

def combine_qsm(in_folder, out_folder):
    """
    1. find all contrasts input folder
    2. call combine_acquisitions for all contrasts
        a. append all readout combined images (accending TE) for phase and magnitude
        b. save multi-TE data in output folder
    """
    files = [f for f in os.listdir(in_folder) if ("e1.nii" in f)]
    for f in files:
        combine_acquisitions(f, in_folder, out_folder)
        print("acq combined %s" %os.path.join(in_folder,f))

def combine_complex(f, in_folder, out_folder):
    # strip ".nii"  from filename
    filename = f[:-4]

    # load data
    mag = nib.load(os.path.join(in_folder,f))
    mag_data = mag.get_fdata().astype(np.float32)

    pha = nib.load(os.path.join(in_folder,filename+"_ph.nii"))
    pha_data = pha.get_fdata().astype(np.float32)

    # create complex dataset
    comp = mag_data * np.exp(2.0j*pha_data * np.pi/4096)

    # phase correction
    pha_diff = np.angle(comp[...,1:]*np.conj(comp[...,:-1]))
    sigma = np.zeros(len(pha_diff.shape))
    sigma[:3] = 2
    corr_real = gaussian_filter(np.cos(pha_diff), sigma=sigma)
    corr_imag = gaussian_filter(np.sin(pha_diff), sigma=sigma)
    corr = np.exp(1.0j*np.arctan2(corr_imag, corr_real))
    comp[...,1:] *= np.conj(corr)

    # average over readout direction
    comp= np.sum(comp,axis=-1)

    #save data
    mag_im = nib.Nifti1Image(abs(comp), mag.affine,mag.header)
    nib.save(mag_im, os.path.join(out_folder, f))

    pha_im = nib.Nifti1Image(np.angle(comp), mag.affine,mag.header)
    nib.save(pha_im, os.path.join(out_folder, filename+"_ph.nii"))

    # copy json files
    os.system("cp %s %s" % (os.path.join(in_folder,filename+".json"), os.path.join(out_folder,filename+".json")))

def combine_acquisitions(f, in_folder, out_folder):
    # create empty lists for magnitude and phase data
    comb_mag = []
    comb_pha = []

    # append data to lists for each TE
    for i in range(1,5):
        filename = f[:-7]
        mag = nib.load(os.path.join(in_folder,filename+"_e%s.nii" %i))
        comb_mag.append(mag.get_fdata().astype(np.float32))
        pha = nib.load(os.path.join(in_folder,filename+"_e%s_ph.nii" %i))
        comb_pha.append(pha.get_fdata().astype(np.float32))

    # save data; note that TE is now the 1st dimension and need to be the last -->moveaxis
    mag_im = nib.Nifti1Image(np.moveaxis(np.array(comb_mag),0,-1), mag.affine,mag.header)
    nib.save(mag_im, os.path.join(out_folder, filename+"_mag_roc.nii"))

    pha_im = nib.Nifti1Image(np.moveaxis(np.array(comb_pha),0,-1), mag.affine,mag.header)
    nib.save(pha_im, os.path.join(out_folder, filename+"_pha_roc.nii"))

    ## todo for Monica: create proper json for this nifti and save it in out_folder

def create_folder(path):
    if not os.path.isdir(path):
        os.mkdir(path)

def combine(path):
    sites    = ["DZNE"]
    subjects = ["subj-01", "subj-02", "subj-03"]
    sessions = ["ptx", "cp"]

    for site in sites:
        for subject in subjects:
            for session in sessions:
                mpm_out = os.path.join(path,"derivatives",site,subject,session,"mpm","ROCombine")
                qsm_out = os.path.join(path,"derivatives",site,subject,session,"qsm","ROCombine")
                create_folder(mpm_out)
                create_folder(qsm_out)
                combine_mpm(os.path.join(path,site,subject,session,"mpm"),mpm_out)
                combine_qsm(mpm_out,qsm_out)

combine("/Users/voelzkey/Desktop/Data/QSMData/2209_talk")