import os
import sys
import argparse


def mpm_ptx(path,site,subject):
    filename_batch = os.path.join(path,"derivatives",site,subject,"ptx","mpm","spm_batch.m")
    input_folder   = os.path.join(path,"derivatives",site,subject,"ptx","mpm","ROcombine")
    output_folder  = os.path.join(path,"derivatives",site,subject,"ptx","mpm","maps")

    f = open(filename_batch,"w")

    # begin file
    f.write("%---------------------------------------\n")
    f.write("% This is an automatically generated batchfile\n")
    f.write("% author: voelzkey\n")
    f.write("% date: xxx\n")
    f.write("%---------------------------------------\n\n")

    f.write("matlabbatch{1}.spm.tools.hmri.create_mpm.subj.output.outdir = {'%s'};\n" %output_folder)
    f.write("matlabbatch{1}.spm.tools.hmri.create_mpm.subj.b1_type.no_B1_correction = 'noB1';\n")
    f.write("matlabbatch{1}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_us = '-';\n")

    # input data MT
    f.write("matlabbatch{1}.spm.tools.hmri.create_mpm.subj.raw_mpm.MT = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_ptx_MTwUP_e1.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_ptx_MTwUP_e2.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_ptx_MTwUP_e3.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_ptx_MTwUP_e4.nii" %(site,subject)))
    f.write("                                                            };\n")

    # input data PD
    f.write("matlabbatch{1}.spm.tools.hmri.create_mpm.subj.raw_mpm.PD = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_ptx_PD_e1.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_ptx_PD_e2.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_ptx_PD_e3.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_ptx_PD_e4.nii" %(site,subject)))
    f.write("                                                            };\n")

    # input data T1
    f.write("matlabbatch{1}.spm.tools.hmri.create_mpm.subj.raw_mpm.T1 = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_ptx_T1_e1.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_ptx_T1_e2.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_ptx_T1_e3.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_ptx_T1_e4.nii" %(site,subject)))
    f.write("                                                            };\n")

    # Disable popups and close file
    f.write("matlabbatch{1}.spm.tools.hmri.create_mpm.subj.popup = false;\n")
    f.close()

    call_batch(filename_batch)

def mpm_cp(path,site,subject):
    filename_batch = os.path.join(path,"derivatives",site,subject,"cp","mpm","spm_batch.m")
    input_folder   = os.path.join(path,"derivatives",site,subject,"cp","mpm","ROcombine")
    output_folder  = os.path.join(path,"derivatives",site,subject,"cp","mpm","maps")
    b1_folder      = os.path.join(path,site,subject,"cp","fmap")

    f = open(filename_batch,"w")

    # begin file
    f.write("%---------------------------------------\n")
    f.write("% This is an automatically generated batchfile\n")
    f.write("% author: voelzkey\n")
    f.write("% date: xxx\n")
    f.write("%---------------------------------------\n\n")

    f.write("matlabbatch{1}.spm.tools.hmri.create_mpm.subj.output.outdir = {'%s'};\n" %output_folder)
    f.write("matlabbatch{1}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_us = '-';\n")

    f.write("matlabbatch{1}.spm.tools.hmri.create_mpm.subj.b1_type.pre_processed_B1.b1input = {\n")
    f.write("                                                                                  '%s,1'\n"%os.path.join(b1_folder,"%s_%s_cp_fmap_b1_con.nii" %(site,subject)))
    f.write("                                                                                  '%s,1'\n"%os.path.join(b1_folder,"%s_%s_cp_fmap_b1.nii" %(site,subject)))
    f.write("                                                                                  };\n")
    f.write("matlabbatch{1}.spm.tools.hmri.create_mpm.subj.b1_type.pre_processed_B1.scafac = 0.1;\n")


    # input data MT
    f.write("matlabbatch{1}.spm.tools.hmri.create_mpm.subj.raw_mpm.MT = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_cp_MTwCP_e1.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_cp_MTwCP_e2.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_cp_MTwCP_e3.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_cp_MTwCP_e4.nii" %(site,subject)))
    f.write("                                                            };\n")

    # input data PD
    f.write("matlabbatch{1}.spm.tools.hmri.create_mpm.subj.raw_mpm.PD = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_cp_PD_e1.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_cp_PD_e2.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_cp_PD_e3.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_cp_PD_e4.nii" %(site,subject)))
    f.write("                                                            };\n")

    # input data T1
    f.write("matlabbatch{1}.spm.tools.hmri.create_mpm.subj.raw_mpm.T1 = {\n")
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_cp_T1_e1.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_cp_T1_e2.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_cp_T1_e3.nii" %(site,subject)))
    f.write("                                                            '%s'\n"%os.path.join(input_folder,"%s_%s_cp_T1_e4.nii" %(site,subject)))
    f.write("                                                            };\n")

    # Disable popups and close file
    f.write("matlabbatch{1}.spm.tools.hmri.create_mpm.subj.popup = false;\n")
    f.close()

    call_batch(filename_batch)

def call_batch(filename):
    import matlab.engine as mat # move matlab import to here if matlab is not installed only this function fails

    eng=mat.start_matlab()
    eng.call_batch(filename,nargout=0)
    eng.quit()

def main():
    parser = argparse.ArgumentParser(
    description='Processing pipeline for MPM/QSM data in SCAIFIELD piloting.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--path',   help='path to MPM/QSM data, folder, where the BIDS structure originates ', required=True)
    parser.add_argument('--site',   help='sites that are supposed to be analysed, default: a for all', required=False, default = 'a', type=list)
    parser.add_argument('--subj',   help='subjects that are supposed to be analysed, default: a for all', required=False, default = 'a', type=list)
    
    args = parser.parse_args()
    path = args.path
    sites = args.site
    subjects = args.subj

    print(sites)
    if sites[0] == "a":
        sites = [ f.name for f in os.scandir(os.path.join(path,"derivatives")) if f.is_dir() ]
    
    print(sites)
    
    for site in sites:
        print(site)
        if subjects[0] == 'a':
            for s in [ f.name for f in os.scandir(os.path.join(path,"derivatives",site)) if f.is_dir() ]:
                mpm_ptx(path,site,s)
                mpm_cp(path,site,s)
        else:
            for s in subjects:
                mpm_ptx(path,site,s)
                mpm_cp(path,site,s)

if __name__ == '__main__':
    sys.exit(main())
