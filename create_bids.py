import os
import numpy as np

"""
This script add one subject (ptx and cp) to the BIDS structure and creates a this structure, if it does not exists

Note that it is based on the names from Bonn XNAT and it might need changes for other sites / data storage methods 
"""
path_bids  = "/Users/voelzkey/Desktop/Data/QSMData/2209_talk"  # folder where you want to build your BIDS structure
path_dcm   = "scans3"                                          # folder with the dicom scans (download from xnat and extract)

site = "DZNE"    # site in BIDS (keep DZNE)
subj = "subj-03" # name of the subject... there is a check if subject already exists

"""
all changes needed are above
"""
# ugly hardcoded stuff
path_site = os.path.join(path_bids,site); 
path_subj = os.path.join(path_site,subj); 
path_der  = os.path.join(path_bids,"derivatives"); 

path_ptx = os.path.join(path_subj,"ptx")
path_magn_ptx = os.path.join(path_ptx,"magn")
path_fmap_ptx = os.path.join(path_ptx,"fmap")
path_mpm_ptx  = os.path.join(path_ptx,"mpm")

path_cp = os.path.join(path_subj,"cp")
path_magn_cp = os.path.join(path_cp,"magn")
path_fmap_cp = os.path.join(path_cp,"fmap")
path_mpm_cp  = os.path.join(path_cp,"mpm")

path_der_site = os.path.join(path_der,site)
path_der_subj = os.path.join(path_der_site,subj)

path_der_ptx = os.path.join(path_der_subj,"ptx")
path_der_magn_ptx = os.path.join(path_der_ptx,"magn")
path_der_fmap_ptx = os.path.join(path_der_ptx,"fmap")
path_der_mpm_ptx  = os.path.join(path_der_ptx,"mpm")
path_der_qsm_ptx  = os.path.join(path_der_ptx,"qsm")

path_der_cp = os.path.join(path_der_subj,"cp")
path_der_magn_cp = os.path.join(path_der_cp,"magn")
path_der_fmap_cp = os.path.join(path_der_cp,"fmap")
path_der_mpm_cp  = os.path.join(path_der_cp,"mpm")
path_der_qsm_cp  = os.path.join(path_der_cp,"qsm")

def create_folder(path):
    if not os.path.isdir(path):
        os.mkdir(path)

def check_and_create_folder(path):
    if os.path.isdir(path):
        raise ValueError("%s alrady exists" %path)
    else:
        os.mkdir(path)


def mkdir():
    create_folder(path_site)
    check_and_create_folder(path_subj) # check if subject already exists and throw error to not overwrite data
    create_folder(path_der)

    create_folder(path_ptx)
    create_folder(path_magn_ptx)
    create_folder(path_fmap_ptx)
    create_folder(path_mpm_ptx)

    create_folder(path_cp)
    create_folder(path_magn_cp)
    create_folder(path_fmap_cp)
    create_folder(path_mpm_cp)

    create_folder(path_der_site)
    create_folder(path_der_subj)

    create_folder(path_der_ptx)
    create_folder(path_der_magn_ptx)
    create_folder(path_der_fmap_ptx)
    create_folder(path_der_mpm_ptx)
    create_folder(path_der_qsm_ptx)

    create_folder(path_der_cp)
    create_folder(path_der_magn_cp)
    create_folder(path_der_fmap_cp)
    create_folder(path_der_mpm_cp)
    create_folder(path_der_qsm_cp)

def parse_raw_data():
    prefix_ptx = "%s_%s_ptx_" %(site,subj)
    prefix_cp  = "%s_%s_cp_" %(site,subj)

    dcm_files = sorted(os.listdir(path_dcm))
    for f in dcm_files:
        if "MTwUP_" in f:
            name = "%sMTwUP" %(prefix_ptx)
            if "E00" in f:
                name+= "_e1"
            folder = path_mpm_ptx
        elif "_ptxmpm_MTwCP" in f:
            name = "%sMTwCP" %(prefix_ptx)
            if "E00" in f:
                name+= "_e1"
            folder = path_mpm_ptx
        elif "_mpm_MTwCP" in f:
            name = "%sMTwCP" %(prefix_cp)
            folder = path_mpm_cp
            if "E00" in f:
                name+= "_e1"
        elif "_ptxmpm_PD" in f:
            name = "%sPD" %(prefix_ptx)
            folder = path_mpm_ptx
            if "E00" in f:
                name+= "_e1"
        elif "_mpm_PD" in f:
            name = "%sPD" %(prefix_cp)
            folder = path_mpm_cp
            if "E00" in f:
                name+= "_e1"
        elif "_ptxmpm_T1" in f:
            name = "%sT1" %(prefix_ptx)
            folder = path_mpm_ptx
            if "E00" in f:
                name+= "_e1"
        elif "_mpm_T1" in f:
            name = "%sT1" %(prefix_cp)
            folder = path_mpm_cp
            if "E00" in f:
                name+= "_e1"
        elif "_rel_B1" in f:
            name = "%sfmap_B1SC" %(prefix_ptx)
            folder = path_fmap_ptx
        elif "B1Pha" in f:
            name = "%sfmap_B1SC" %(prefix_ptx)
            folder = path_fmap_ptx
        elif "B1Com" in f:
            name = "%sfmap_B1" %(prefix_cp)
            folder = path_fmap_cp
        else:
            continue
        os.system("dcm2niix -f %s -o %s %s" %(name, folder, os.path.join(path_dcm,f)))
        #print(name)

def ROcombine_MPM():
    pass

def ROcombine_QSM():
    pass

mkdir()
parse_raw_data()

