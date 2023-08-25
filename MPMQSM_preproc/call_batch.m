function call_batch(path_to_batch)

    %%% SITE-DEPENDENT VARIABLES - TO BE CHANGED %%%
    hmri_path = '/home/marcantf/software/hMRI-toolbox-0.3.0';
    SPM_path  = '/home/marcantf/software/spm12';

    addpath(hmri_path);
    addpath(genpath(SPM_path));
    disp("path addded")

    run(path_to_batch)
    disp(matlabbatch)
    spm_jobman('run',matlabbatch)
end
