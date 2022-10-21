import os
import numpy as np
import nibabel as nib
import json
import sys
import argparse
from sim_FA import Pulse_from_ini

from scipy.ndimage import gaussian_filter
import shutil

"""
This scripts creates readout combined images from the raw MPM/QSM acquisitions for ALL sites/subjects/sessions/acquisition
It saves these data in the format that is needed for MPM and QSM in derivatives/subject/site/session/mpm and derivatives/subject/site/session/mpm, respectively

change path to BIDS in the very last line of this file
"""

def combine_mpm(in_folder, out_folder, raw_folder):
    """
    1. find all pairs of magnitide/phase images in input folder
    2. call combine_complex for all pairs
        a. combine readout directions
        b. save readout combined images to output folder
    """
    files = [f for f in os.listdir(in_folder) if ((".nii" in f) and not ("ph.nii" in f))]
    for f in files:
        combine_complex(f, in_folder, out_folder, raw_folder)
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

def combine_complex(f, in_folder, out_folder, raw_folder):
    
    # strip ".nii" or ".nii.gz" from filename
    filename = f[:-4]
    if ".nii.gz" in f:
        filename = filename[:-3]
    
    # load data
    mag = nib.load(os.path.join(in_folder,f))
    mag_data = mag.get_fdata().astype(np.float32)

    
    if ".nii.gz" in f:
        pha = nib.load(os.path.join(in_folder,filename+"_ph.nii.gz"))
    else:
        pha = nib.load(os.path.join(in_folder,filename+"_ph.nii"))
    pha_data = pha.get_fdata().astype(np.float32)

    if not "_den" in filename:
        pha_data *= np.pi/4096
    # create complex dataset
    comp = mag_data * np.exp(2.0j*pha_data)

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

    if "_den" in filename:
        filename = filename[:-4]
        json_folder = raw_folder
    else:
        json_folder =in_folder

    #save data
    mag_im = nib.Nifti1Image(abs(comp), mag.affine,mag.header)
    nib.save(mag_im, os.path.join(out_folder, filename + ".nii"))

    pha_im = nib.Nifti1Image(np.angle(comp), mag.affine,mag.header)
    nib.save(pha_im, os.path.join(out_folder, filename+"_ph.nii"))

    # correct timings (ms/s conversion) in json files and save in mpm-output folder
    with open(os.path.join(json_folder,filename+".json"),"r+") as f:
        param = json.load(f)
        param["EchoTime"] = 1000*param["EchoTime"] 
        param["RepetitionTime"] = 1000*param["RepetitionTime"] 

    with open(os.path.join(out_folder,filename+".json"),"w") as f:
        json.dump(param,f,indent=4)

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

def denoise(in_folder, out_folder,prelim):
    create_folder(out_folder)
    create_folder(prelim)

    """files = [f for f in os.listdir(in_folder) if ((".nii" in f) and not ("ph.nii" in f))]
    is_ptx = "_ptx_" in files[0]
    if is_ptx:
        files = [f for f in files if "MTwCP" not in f]
    
    command_magn = "fslmerge -t %s/magn.nii" %prelim
    for f in files:
        command_magn += " %s" % os.path.join(in_folder,f)
    os.system(command_magn)
    
    command_pha = "fslmerge -t %s/pha.nii" %prelim
    for f in files:
        command_pha += " %s" % os.path.join(in_folder,f[:-4]+"_ph.nii")
    os.system(command_pha)

    os.system("mrcalc %s/pha.nii.gz 0.000766990393943 -mult %s/tmp_pha.nii.gz" %(prelim,prelim))

    os.system("mrcalc %s/magn.nii.gz %s/tmp_pha.nii.gz -polar %s/tmp_comp.mif --force" %(prelim,prelim,prelim))

    os.system("dwidenoise %s/tmp_comp.mif %s/out.mif" %(prelim,prelim))

    os.system("mrcalc %s/out.mif -abs %s/out_magn.nii --force" %(prelim,prelim))
    os.system("mrcalc %s/out.mif -phase %s/out_phas.nii --force" %(prelim,prelim))
    
    for i in range(len(files)):
        os.system("fslroi %s/out_magn.nii %s %i 2" %(prelim,os.path.join(out_folder,files[i][:-4]+"_den.nii"), 2*i))
        os.system("fslroi %s/out_phas.nii %s %i 2" %(prelim,os.path.join(out_folder,files[i][:-4]+"_den_ph.nii"), 2*i))

    shutil.rmtree(prelim, ignore_errors=True)"""
    #os.rmdir(prelim)

def calc_FAmap(ini_file, b1_in_file, b1_out_file):

        Pulse = Pulse_from_ini(ini_file,b1_in_file)

        X,Y,Z,I = Pulse.B1.shape

        FA = np.zeros((X,Y,Z))
        PH = np.zeros((X,Y,Z))

        for x in range(X):
            for y in range(Y):
                for z in range(Z):
                    print(x,y,z)
                    FA[x,y,z],PH[x,y,z] = Pulse.calc_FA(x,y,z,0)

        FA *= 100 / Pulse.nomFA

        mag_im = nib.Nifti1Image(FA, Pulse.affine,Pulse.header)
        nib.save(mag_im, b1_out_file)

def main():
    parser = argparse.ArgumentParser(
    description='Processing pipeline for MPM/QSM data in SCAIFIELD piloting.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--path', help='path to MPM/QSM data, folder, where the BIDS structure originates ', required=True)
    parser.add_argument('--den',  help='boolean, shall the raw data be denoised? default = true ', required=False, default=True)

    args = parser.parse_args()
    path = args.path
    den  = args.den

    sites = [ f.name for f in os.scandir(os.path.join(path,"derivatives")) if f.is_dir() ]
    
    print(sites)
    
    for site in sites:
        print(site)
        for subject in sorted([ f.name for f in os.scandir(os.path.join(path,"derivatives",site)) if f.is_dir()]):
            for session in [f.name for f in os.scandir(os.path.join(path,"derivatives",site,subject)) if f.is_dir()]:
                if "ptx" in session:
                    b1_in    = os.path.join(path,site,subject,session,"fmap","%s_%s_ptx_fmap_B1SC.nii" %(site,subject))
                    b1_out   = os.path.join(path,"derivatives",site,subject,session,"fmap","%s_%s_ptx_fmap_B1_sim.nii" %(site,subject))
                    ini_file = "pTXNormal.ini"

                    #calc_FAmap(ini_file, b1_in, b1_out)

                raw_folder = os.path.join(path,site,subject,session,"mpm")
                if den:
                    mpm_in  = os.path.join(path,"derivatives",site,subject,session,"mpm","denoised")
                    prelim  = os.path.join(path,"derivatives",site,subject,session,"mpm","prelim")
                    denoise(raw_folder,mpm_in,prelim)
                else:

                    mpm_in = raw_folder                
                mpm_out = os.path.join(path,"derivatives",site,subject,session,"mpm","ROCombine")
                qsm_out = os.path.join(path,"derivatives",site,subject,session,"qsm","ROCombine")
                create_folder(mpm_out)
                create_folder(qsm_out)
                combine_mpm(mpm_in,mpm_out,raw_folder)
                combine_qsm(mpm_out,qsm_out)

if __name__ == '__main__':
    sys.exit(main())