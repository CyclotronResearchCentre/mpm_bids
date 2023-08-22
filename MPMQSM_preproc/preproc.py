import os
import numpy as np
import nibabel as nib
import json
import sys
import argparse

from scipy.ndimage import gaussian_filter
import shutil

"""
This scripts creates readout combined images from the raw MPM/QSM acquisitions
It saves these data in the format that is needed for MPM and QSM in derivatives/subject/site/session/mpm and derivatives/subject/site/session/mpm, respectively
"""

def combine_mpm(in_folder, out_folder, raw_folder):
    """
    1. find all pairs of magnitide/phase images in input folder
    2. call combine_complex for all pairs
        a. combine readout directions
        b. save readout combined images to output folder
    """
    files = [f for f in os.listdir(in_folder) if ((".nii" in f) and not ("ph.nii" in f) and not ("._" in f))]
    for f in files:
        combine_complex(f, in_folder, out_folder, raw_folder)
        print("RO combined %s" %os.path.join(in_folder,f))

def combine_qsm(in_folder, out_folder):
    """
    1. find all contrasts input folder
    2. call combine_acquisitions for all contrasts
        a. append all readout combined images (ascending TE) for phase and magnitude
        b. save multi-TE data in output folder
    """
    files = [f for f in os.listdir(in_folder) if ("e1.nii" in f)]
    files = [f for f in files if ("._" not in f)]
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
    comp = mag_data * np.exp(1.0j*pha_data)

    # outcomment phase complex averaging, as it is not performed anymore
    """# phase correction
    pha_diff = np.angle(comp[...,1:]*np.conj(comp[...,:-1]))
    sigma = np.zeros(len(pha_diff.shape))
    sigma[:3] = 2
    corr_real = gaussian_filter(np.cos(pha_diff), sigma=sigma)
    corr_imag = gaussian_filter(np.sin(pha_diff), sigma=sigma)
    corr = np.exp(1.0j*np.arctan2(corr_imag, corr_real))
    comp[...,1:] *= np.conj(corr)

    # average over readout direction
    comp= np.sum(comp,axis=-1)"""

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

    """
    # correct timings (ms/s conversion) in json files and save in mpm-output folder
    with open(os.path.join(json_folder,filename+".json"),"r+") as f:
        param = json.load(f)
        param["EchoTime"] = 1000*param["EchoTime"] 
        param["RepetitionTime"] = 1000*param["RepetitionTime"] 

    with open(os.path.join(out_folder,filename+".json"),"w") as f:
        json.dump(param,f,indent=4)
    """

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
    print(path)
    if not os.path.isdir(path):
        os.mkdir(path)

def denoise(in_folder, out_folder,prelim):
    print(out_folder, prelim)
    create_folder(out_folder)
    create_folder(prelim)

    files = [f for f in os.listdir(in_folder) if ((".nii" in f) and not ("ph.nii" in f))]
    print(files)
    
    command_magn = "fslmerge -t %s/magn.nii" %prelim
    for f in files:
        command_magn += " %s" % os.path.join(in_folder,f)
    print(command_magn)
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
        os.system("fslroi %s/out_magn.nii %s %i 1" %(prelim,os.path.join(out_folder,files[i][:-4]+"_den.nii"), i))
        os.system("fslroi %s/out_phas.nii %s %i 1" %(prelim,os.path.join(out_folder,files[i][:-4]+"_den_ph.nii"), i))

    shutil.rmtree(prelim, ignore_errors=True)

def combine_acquisitions_moco(f, in_folder, out_folder):
    # create empty lists for magnitude and phase data
    comb_mag = []
    comb_pha = []

    # append data to lists for each TE
    for i in range(1,5):
        filename = f[:-7]
        real = nib.load(os.path.join(in_folder,filename+"_e%s_real_moco.nii.gz" %i))
        mag  = real.get_fdata().astype(np.complex64)
        imag = nib.load(os.path.join(in_folder,filename+"_e%s_imag_moco.nii.gz" %i))
        mag += 1j*imag.get_fdata().astype(np.float32)
        
        comb_mag.append(abs(mag))
        comb_pha.append(np.angle(mag))

    # save data; note that TE is now the 0st dimension and need to be the last -->moveaxis
    mag_im = nib.Nifti1Image(np.moveaxis(np.array(comb_mag),0,-1), real.affine,real.header)
    nib.save(mag_im, os.path.join(out_folder, filename+"_mag_roc.nii"))

    pha_im = nib.Nifti1Image(np.moveaxis(np.array(comb_pha),0,-1), real.affine,real.header)
    nib.save(pha_im, os.path.join(out_folder, filename+"_pha_roc.nii"))

def moco(mpm_path,moco_path,qsm_path):
    create_folder(moco_path)
    files = os.listdir(mpm_path)
    print(files)
    print(moco_path, qsm_path)
    first_im = sorted([f for f in files if "e1.nii" in f] )
    if len(first_im) == 0:
        first_im = sorted([f for f in files if "e1_den.nii" in f] )


    in_file  = os.path.join(mpm_path, first_im[1])
    ref_file = os.path.join(mpm_path, first_im[2])
    out_file = os.path.join(moco_path, first_im[1][:-4])

    commandPD = "flirt -in %s -out %s -ref %s -omat PD.mat" %(in_file, out_file, ref_file)

    in_file  = os.path.join(mpm_path, first_im[0])
    ref_file = os.path.join(mpm_path, first_im[2])
    out_file = os.path.join(moco_path, first_im[0][:-4])
    commandMT = "flirt -in %s -out %s -ref %s -omat MT.mat" %(in_file, out_file, ref_file)

    os.system(commandMT)
    os.system(commandPD)

    PD_files = [f for f in files if "_pd_" in f]
    PD_files = sorted([f[:-4] for f in PD_files if ((".nii" in f) and ("ph.nii" not in f))])

    for f in PD_files:
        im = nib.load(os.path.join(mpm_path,f+".nii"))
        mag = im.get_fdata()
        pha = nib.load(os.path.join(mpm_path,f+"_ph.nii")).get_fdata()
        
        com = mag*np.exp(1j*pha)
        nib.save(nib.Nifti1Image(com.real, affine = im.affine, header=im.header), os.path.join(moco_path,f+"_real.nii"))
        nib.save(nib.Nifti1Image(com.imag, affine = im.affine, header=im.header), os.path.join(moco_path,f+"_imag.nii"))
    
    MT_files = [f for f in files if "_mt" in f]
    MT_files = sorted([f[:-4] for f in MT_files if ((".nii" in f) and ("ph.nii" not in f))])

    for f in MT_files:
        im = nib.load(os.path.join(mpm_path,f+".nii"))
        mag = im.get_fdata()
        pha = nib.load(os.path.join(mpm_path,f+"_ph.nii")).get_fdata()
        
        com = mag*np.exp(1j*pha)
        nib.save(nib.Nifti1Image(com.real, affine = im.affine, header=im.header), os.path.join(moco_path,f+"_real.nii"))
        nib.save(nib.Nifti1Image(com.imag, affine = im.affine, header=im.header), os.path.join(moco_path,f+"_imag.nii"))

    ref_file = os.path.join(mpm_path, first_im[2])

    for f in PD_files:
        in_file  = os.path.join(moco_path, f+"_real")
        out_file = os.path.join(moco_path, f+"_real_moco")
        command  = "flirt -ref %s -in %s -out %s -applyxfm -init PD.mat" % (ref_file, in_file, out_file)
        os.system(command)
        
        in_file  = os.path.join(moco_path, f+"_imag")
        out_file = os.path.join(moco_path, f+"_imag_moco")
        command  = "flirt -ref %s -in %s -out %s -applyxfm -init PD.mat" % (ref_file, in_file, out_file)
        os.system(command)
        
    for f in MT_files:
        in_file  = os.path.join(moco_path, f+"_real")
        out_file = os.path.join(moco_path, f+"_real_moco")
        command  = "flirt -ref %s -in %s -out %s -applyxfm -init MT.mat" % (ref_file, in_file, out_file)
        os.system(command)
        
        in_file  = os.path.join(moco_path, f+"_imag")
        out_file = os.path.join(moco_path, f+"_imag_moco")
        command  = "flirt -ref %s -in %s -out %s -applyxfm -init MT.mat" % (ref_file, in_file, out_file)
        os.system(command)

    combine_acquisitions_moco(first_im[0], moco_path, qsm_path)
    combine_acquisitions_moco(first_im[1], moco_path, qsm_path)


def main():
    parser = argparse.ArgumentParser(
    description='Processing pipeline for MPM/QSM data in SCAIFIELD piloting.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--path', help='path to MPM/QSM data, folder, where the BIDS structure originates ', required=True)
    parser.add_argument('--den',  help='boolean, shall the raw data be denoised? default = true ', required=False, default=True)
    parser.add_argument('--site', help='exact name of site (eg DZNE)', required=True)
    parser.add_argument('--sub',  help='exact name of subject (eg subj-01)', required=True)
    parser.add_argument('--ses',  help='exact name of session (eg ses-01)',  required=True)
    parser.add_argument('--ptx',  help='create log file? default true',   default=True, action='store_true')

    args = parser.parse_args()
    path = args.path
    den  = args.den
    site = args.site
    sub  = args.sub
    ses  = args.ses
    ptx  = args.ptx
    
    raw_folder = os.path.join(path,site,sub,ses,"mpm")
    create_folder(os.path.join(path,"derivatives"))
    create_folder(os.path.join(path,"derivatives",site))
    create_folder(os.path.join(path,"derivatives",site,sub))
    create_folder(os.path.join(path,"derivatives",site,sub,ses))
    create_folder(os.path.join(path,"derivatives",site,sub,ses,"mpm"))
    create_folder(os.path.join(path,"derivatives",site,sub,ses,"qsm"))
    if den:
        mpm_in  = os.path.join(path,"derivatives",site,sub,ses,"mpm","denoised")
        prelim  = os.path.join(path,"derivatives",site,sub,ses,"mpm","prelim")
        #denoise(raw_folder,mpm_in,prelim)
    else:
        mpm_in = raw_folder    
         
    mpm_out = os.path.join(path,"derivatives",site,sub,ses,"mpm","ROCombine")
    qsm_out = os.path.join(path,"derivatives",site,sub,ses,"qsm","ROCombine")
    create_folder(mpm_out)
    create_folder(qsm_out)
    combine_mpm(mpm_in,mpm_out,raw_folder)
    combine_qsm(mpm_out,qsm_out)
    
    path_moco = os.path.join(path,"derivatives",site,sub,ses,"qsm","moco")
    moco(mpm_out,path_moco,qsm_out)
    
    mpm(path, site, sub, ses, ptx)

def call_batch(filename):
    import matlab.engine as mat # move matlab import to here if matlab is not installed only this function fails

    eng=mat.start_matlab()
    eng.addpath(os.path.dirname(filename))
    eng.addpath(os.path.dirname(__file__))
    print(os.path.dirname(filename))
    eng.call_batch(filename,nargout=0)
    eng.quit()

def mpm_ptx(path,site,subject,session):
    
    filename_batch = os.path.join(path,"derivatives",site,subject,session,"mpm","spm_batch.m")
    input_folder   = os.path.join(path,"derivatives",site,subject,session,"mpm","ROcombine")
    output_folder  = os.path.join(path,"derivatives",site,subject,session,"mpm","maps")
    b1_folder      = os.path.join(path,"derivatives",site,subject,session,"fmap")
    b1_raw         = os.path.join(path,site,subject,session,"fmap")

    f = open(filename_batch,"w")

    # begin file
    f.write("%---------------------------------------\n")
    f.write("% This is an automatically generated batchfile\n")
    f.write("% author: voelzkey\n")
    f.write("% date: xxx\n")
    f.write("%---------------------------------------\n\n")

    f.write("matlabbatch{1}.spm.tools.hmri.hmri_config.hmri_setdef.customised = {'%s'};\n" % os.path.join(os.path.dirname(os.path.abspath(__file__)),"hmri_defaults_SS_scaifield.m"))

    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.output.outdir = {'%s'};\n" %output_folder)
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.pre_processed_B1.b1input = {\n")
    f.write("                                                                                  '%s,1'\n"%os.path.join(b1_raw,"%s_%s_%s_fmap-b1-con.nii" %(site,subject,session)))
    f.write("                                                                                  '%s,1'\n"%os.path.join(b1_raw,"%s_%s_%s_fmap-b1.nii" %(site,subject,session)))
    f.write("                                                                                  };\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.pre_processed_B1.scafac = .1;\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_us = '-';\n")

    # input data MT
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.MT = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_mt_e1.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_mt_e2.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_mt_e3.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_mt_e4.nii" %(site,subject,session)))
    f.write("                                                            };\n")

    # input data PD
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.PD = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_pd_e1.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_pd_e2.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_pd_e3.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_pd_e4.nii" %(site,subject,session)))
    f.write("                                                            };\n")

    # input data T1
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.T1 = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_t1_e1.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_t1_e2.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_t1_e3.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_t1_e4.nii" %(site,subject,session)))
    f.write("                                                            };\n")

    # Disable popups and close file
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.popup = false;\n")
    f.close()

    print(filename_batch)
    call_batch(filename_batch)

def mpm_cp(path,site,subject,session):
    filename_batch = os.path.join(path,"derivatives",site,subject,session,"mpm","spm_batch.m")
    input_folder   = os.path.join(path,"derivatives",site,subject,session,"mpm","ROcombine")
    output_folder  = os.path.join(path,"derivatives",site,subject,session,"mpm","maps")
    b1_folder      = os.path.join(path,site,subject,session,"fmap")

    f = open(filename_batch,"w")

    # begin file
    f.write("%---------------------------------------\n")
    f.write("% This is an automatically generated batchfile\n")
    f.write("% author: voelzkey\n")
    f.write("% date: xxx\n")
    f.write("%---------------------------------------\n\n")

    f.write("matlabbatch{1}.spm.tools.hmri.hmri_config.hmri_setdef.customised = {'%s'};\n" % os.path.join(os.path.dirname(os.path.abspath(__file__)),"hmri_defaults_SS_scaifield.m"))

    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.output.outdir = {'%s'};\n" %output_folder)
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_us = '-';\n")

    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.pre_processed_B1.b1input = {\n")
    f.write("                                                                                  '%s,1'\n"%os.path.join(b1_folder,"%s_%s_%s_fmap_b1_con.nii" %(site,subject,session)))
    f.write("                                                                                  '%s,1'\n"%os.path.join(b1_folder,"%s_%s_%s_fmap_b1.nii" %(site,subject,session)))
    f.write("                                                                                  };\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.pre_processed_B1.scafac = %f;\n" %(335/2400)) #correct for differences in refVol


    # input data MT
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.MT = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_mt_e1.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_mt_e2.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_mt_e3.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_mt_e4.nii" %(site,subject,session)))
    f.write("                                                            };\n")

    # input data PD
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.PD = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_pd_e1.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_pd_e2.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_pd_e3.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_pd_e4.nii" %(site,subject,session)))
    f.write("                                                            };\n")

    # input data T1
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.T1 = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_t1_e1.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_t1_e2.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_t1_e3.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_t1_e4.nii" %(site,subject,session)))
    f.write("                                                            };\n")

    # Disable popups and close file
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.popup = false;\n")
    f.close()

    print(filename_batch)
    call_batch(filename_batch)

def mpm(path, site, subject, session, isptxdata):
    if isptxdata:
        mpm_ptx(path,site,subject,session)
    else:
        mpm_cp(path,site,subject,session)

if __name__ == '__main__':
    sys.exit(main())