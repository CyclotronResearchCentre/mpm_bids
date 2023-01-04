function call_batch(path_to_batch)
    hmri_path = '/Users/voelzkey/Desktop/CodeMatlab/hMRI-toolbox-0.2.4';
    SPM_path  = '/Users/voelzkey/Desktop/CodeMatlab/spm12';

    addpath(hmri_path);
    addpath(genpath(SPM_path));
    disp("path addded")

    run(path_to_batch)
    disp(matlabbatch)
    spm_jobman('run',matlabbatch)
end
