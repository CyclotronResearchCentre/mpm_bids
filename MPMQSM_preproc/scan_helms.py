import os
import numpy as np
import nibabel as nib
import json
import sys
import argparse

from scipy.ndimage import gaussian_filter
import shutil
import configparser

def create_folder(path):
    print(path)
    if not os.path.isdir(path):
        os.mkdir(path)



def main():
    parser = argparse.ArgumentParser(
    description='scan helms',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--path', help='path to MPM/QSM data, folder, where the BIDS structure originates ', required=True)
    parser.add_argument('--den',  help='boolean, shall the raw data be denoised? default = true ', required=False, default=True)
    parser.add_argument('--site', help='exact name of site (eg DZNE)', required=True)
    parser.add_argument('--sub',  help='exact name of subject (eg subj-01)', required=True)
    parser.add_argument('--ses',  help='exact name of session (eg ses-01)',  required=True)
    parser.add_argument('--ptx',  help='create log file? default true', action='store_true')
    parser.add_argument('--name', help='name of the mpm subfolder, default = mpm',default = "mpm",required=False)

    args = parser.parse_args()
    path = args.path
    den  = args.den
    site = args.site
    sub  = args.sub
    ses  = args.ses
    ptx  = args.ptx

    name = args.name
    
    mpm(path, site, sub, ses, ptx, name)

def call_batch(filename):
    import matlab.engine as mat # move matlab import to here if matlab is not installed only this function fails

    eng=mat.start_matlab()
    eng.addpath(os.path.dirname(filename))
    eng.addpath(os.path.dirname(__file__))
    print(os.path.dirname(filename))
    eng.call_batch(filename,nargout=0)
    eng.quit()

def mpm_ptx(path,site,subject,session,name,C,cor):
    if cor == "helms": 
        filename_batch = os.path.join(path,"derivatives",site,subject,session,name,"spm_batch_helms_0p%i.m"%C)
    elif cor == "lipp": 
        filename_batch = os.path.join(path,"derivatives",site,subject,session,name,"spm_batch_%s.m"%C)
    else:
        print("unkown correction method")
        lf
    input_folder   = os.path.join(path,"derivatives",site,subject,session,name,"ROcombine")
    b1_folder      = os.path.join(path,"derivatives",site,subject,session,"fmap")
    b1_raw         = os.path.join(path,site,subject,session,"fmap")

    f = open(filename_batch,"w")

    # begin file
    f.write("%---------------------------------------\n")
    f.write("% This is an automatically generated batchfile\n")
    f.write("% author: voelzkey\n")
    f.write("% date: xxx\n")
    f.write("%---------------------------------------\n\n")

    if cor == "helms": 
        f.write("matlabbatch{1}.spm.tools.hmri.hmri_config.hmri_setdef.customised = {'%s'};\n" % os.path.join(os.path.dirname(os.path.abspath(__file__)),"hmri_defaults/helms_0p%i.m"%C))
        output_folder  = os.path.join(path,"derivatives",site,subject,session,name,"maps_helms_0p%i"%C)
    elif cor == "lipp":
        f.write("matlabbatch{1}.spm.tools.hmri.hmri_config.hmri_setdef.customised = {'%s'};\n" % os.path.join(os.path.dirname(os.path.abspath(__file__)),"hmri_defaults/lipp_%s.m"%C))
        output_folder  = os.path.join(path,"derivatives",site,subject,session,name,"maps_lipp_%s"%C)

    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.output.outdir = {'%s'};\n" %output_folder)
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type.pre_processed_B1.b1input = {\n")
    f.write("                                                                                                '%s,1'\n"%os.path.join(b1_raw,"%s_%s_%s_fmap-afi-con.nii" %(site,subject,session)))
    f.write("                                                                                                '%s,1'\n"%os.path.join(b1_raw,"%s_%s_%s_fmap-afi.nii" %(site,subject,session)))
    f.write("                                                                                  };\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type.pre_processed_B1.scafac = .1;\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type.pre_processed_B1.b1parameters.b1metadata = 'yes';\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type_MT.pre_processed_B1.b1input = {\n")
    f.write("                                                                                                '%s,1'\n"%os.path.join(b1_raw,"%s_%s_%s_fmap-B1SC-con.nii" %(site,subject,session)))
    f.write("                                                                                                '%s,1'\n"%os.path.join(b1_folder,"%s_%s_%s_fmap-B1SC_b1rms.nii" %(site,subject,session)))
    f.write("                                                                                  };\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type_MT.pre_processed_B1.scafac = 100;\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type_MT.pre_processed_B1.b1metadata = 'yes';\n")

    # input data MT
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.MT = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_MTw_e1.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_MTw_e2.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_MTw_e3.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_MTw_e4.nii" %(site,subject,session)))
    f.write("                                                            };\n")

    # input data PD
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.PD = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_PD_e1.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_PD_e2.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_PD_e3.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_PD_e4.nii" %(site,subject,session)))
    f.write("                                                            };\n")

    # input data T1
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.T1 = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_T1_e1.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_T1_e2.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_T1_e3.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_T1_e4.nii" %(site,subject,session)))
    f.write("                                                            };\n")

    # Disable popups and close file
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.popup = false;\n")
    f.close()

    print(filename_batch)
    #call_batch(filename_batch)

def mpm_cp(path,site,subject,session,name,C,cor):
    if cor == "helms": 
        filename_batch = os.path.join(path,"derivatives",site,subject,session,name,"spm_batch_helms_0p%i.m"%C)
    elif cor == "lipp": 
        filename_batch = os.path.join(path,"derivatives",site,subject,session,name,"spm_batch_%s.m"%C)
    else:
        print("unkown correction method")
    input_folder   = os.path.join(path,"derivatives",site,subject,session,name,"ROcombine")
    b1_folder      = os.path.join(path,site,subject,session,"fmap")

    f = open(filename_batch,"w")

    # begin file
    f.write("%---------------------------------------\n")
    f.write("% This is an automatically generated batchfile\n")
    f.write("% author: voelzkey\n")
    f.write("% date: xxx\n")
    f.write("%---------------------------------------\n\n")
    
    if cor == "helms": 
        f.write("matlabbatch{1}.spm.tools.hmri.hmri_config.hmri_setdef.customised = {'%s'};\n" % os.path.join(os.path.dirname(os.path.abspath(__file__)),"hmri_defaults/helms_0p%i.m"%C))
        output_folder  = os.path.join(path,"derivatives",site,subject,session,name,"maps_helms_0p%i"%C)
    elif cor == "lipp":
        f.write("matlabbatch{1}.spm.tools.hmri.hmri_config.hmri_setdef.customised = {'%s'};\n" % os.path.join(os.path.dirname(os.path.abspath(__file__)),"hmri_defaults/lipp_%s.m"%C))
        output_folder  = os.path.join(path,"derivatives",site,subject,session,name,"maps_lipp_%s"%C)

    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.output.outdir = {'%s'};\n" %output_folder)
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_us = '-';\n")

    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.output.outdir = {'%s'};\n" %output_folder)
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type.pre_processed_B1.b1input = {\n")
    f.write("                                                                                  '%s,1'\n"%os.path.join(b1_folder,"%s_%s_%s_fmap-B1SC-con.nii" %(site,subject,session)))
    f.write("                                                                                  '%s,1'\n"%os.path.join(b1_folder,"%s_%s_%s_fmap-B1SC-comb.nii" %(site,subject,session)))
    f.write("                                                                                  };\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type.pre_processed_B1.scafac = .1;\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type.pre_processed_B1.b1parameters.b1metadata = 'yes';\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type_MT.pre_processed_B1.b1input = {\n")
    f.write("                                                                                  '%s,1'\n"%os.path.join(b1_folder,"%s_%s_%s_fmap-B1SC-con.nii" %(site,subject,session)))
    f.write("                                                                                  '%s,1'\n"%os.path.join(b1_folder,"%s_%s_%s_fmap-B1SC-comb.nii" %(site,subject,session)))
    f.write("                                                                                  };\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type_MT.pre_processed_B1.scafac = .1;\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type_MT.pre_processed_B1.b1metadata = 'yes';\n")

    # input data MT
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.MT = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_MTw_e1.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_MTw_e2.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_MTw_e3.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_MTw_e4.nii" %(site,subject,session)))
    f.write("                                                            };\n")

    # input data PD
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.PD = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_PD_e1.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_PD_e2.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_PD_e3.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_PD_e4.nii" %(site,subject,session)))
    f.write("                                                            };\n")

    # input data T1
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.T1 = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_T1_e1.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_T1_e2.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_T1_e3.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_T1_e4.nii" %(site,subject,session)))
    f.write("                                                            };\n")

    # Disable popups and close file
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.popup = false;\n")
    f.close()

    print(filename_batch)
    #call_batch(filename_batch)

def mpm_cp_scaledB1(path,site,subject,session,name,C,cor):
    if cor == "helms": 
        filename_batch = os.path.join(path,"derivatives",site,subject,session,name,"spm_batch_helms_0p%i.m"%C)
    elif cor == "lipp": 
        filename_batch = os.path.join(path,"derivatives",site,subject,session,name,"spm_batch_%s.m"%C)
    else:
        print("unkown correction method")
    input_folder   = os.path.join(path,"derivatives",site,subject,session,name,"ROcombine")
    b1_folder      = os.path.join(path,site,subject,session,"fmap")

    f = open(filename_batch,"w")

    # begin file
    f.write("%---------------------------------------\n")
    f.write("% This is an automatically generated batchfile\n")
    f.write("% author: voelzkey\n")
    f.write("% date: xxx\n")
    f.write("%---------------------------------------\n\n")
    
    if cor == "helms": 
        f.write("matlabbatch{1}.spm.tools.hmri.hmri_config.hmri_setdef.customised = {'%s'};\n" % os.path.join(os.path.dirname(os.path.abspath(__file__)),"hmri_defaults/helms_0p%i.m"%C))
        output_folder  = os.path.join(path,"derivatives",site,subject,session,name,"maps_helms_0p%i_scaledB1"%C)
    elif cor == "lipp":
        f.write("matlabbatch{1}.spm.tools.hmri.hmri_config.hmri_setdef.customised = {'%s'};\n" % os.path.join(os.path.dirname(os.path.abspath(__file__)),"hmri_defaults/lipp_%s.m"%C))
        output_folder  = os.path.join(path,"derivatives",site,subject,session,name,"maps_lipp_%s_scaledB1"%C)

    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.output.outdir = {'%s'};\n" %output_folder)
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_us = '-';\n")

    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.output.outdir = {'%s'};\n" %output_folder)
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type.pre_processed_B1.b1input = {\n")
    f.write("                                                                                  '%s,1'\n"%os.path.join(b1_folder,"%s_%s_%s_fmap-B1SC-con.nii" %(site,subject,session)))
    f.write("                                                                                  '%s,1'\n"%os.path.join(b1_folder,"%s_%s_%s_fmap-B1SC-comb.nii" %(site,subject,session)))
    f.write("                                                                                  };\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type.pre_processed_B1.scafac = .1;\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type.pre_processed_B1.b1parameters.b1metadata = 'yes';\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type_MT.pre_processed_B1.b1input = {\n")
    f.write("                                                                                  '%s,1'\n"%os.path.join(b1_folder,"%s_%s_%s_fmap-B1SC-con.nii" %(site,subject,session)))
    f.write("                                                                                  '%s,1'\n"%os.path.join(b1_folder,"%s_%s_%s_fmap-B1SC-comb.nii" %(site,subject,session)))
    f.write("                                                                                  };\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type_MT.pre_processed_B1.scafac = .132;\n")
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.b1_type.b1_MT.b1_type_MT.pre_processed_B1.b1metadata = 'yes';\n")

    # input data MT
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.MT = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_MTw_e1.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_MTw_e2.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_MTw_e3.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_MTw_e4.nii" %(site,subject,session)))
    f.write("                                                            };\n")

    # input data PD
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.PD = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_PD_e1.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_PD_e2.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_PD_e3.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_PD_e4.nii" %(site,subject,session)))
    f.write("                                                            };\n")

    # input data T1
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.raw_mpm.T1 = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_T1_e1.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_T1_e2.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_T1_e3.nii" %(site,subject,session)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_%s_mpm_T1_e4.nii" %(site,subject,session)))
    f.write("                                                            };\n")

    # Disable popups and close file
    f.write("matlabbatch{2}.spm.tools.hmri.create_mpm.subj.popup = false;\n")
    f.close()

    print(filename_batch)
    #call_batch(filename_batch)

def mpm(path, site, subject, session, isptxdata, name):
    cor = "lipp"
    c = ["1p0","1p2","1p4","1p6","1p8","2p0","2p2","2p4","2p6"]
    for C in c:
        if isptxdata:
            print("ptx")
            mpm_ptx(path,site,subject,session,name,C,cor)
        else:
            print("cp")
            mpm_cp_scaledB1(path,site,subject,session,name,C,cor)
    
    cor = "helms"
    c = [0,1,2,3,4,5,6,7]
    for C in c:
        if isptxdata:
            print("ptx")
            mpm_ptx(path,site,subject,session,name,C,cor)
        else:
            print("cp")
            mpm_cp_scaledB1(path,site,subject,session,name,C,cor)

if __name__ == '__main__':
    sys.exit(main())