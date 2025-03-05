function call_batch(path_to_batch)
    [filepath,name,ext] = fileparts(mfilename("fullpath"));
    hmri_path = fullfile(filepath, 'hMRI-toolbox');
    spm_path = fullfile(filepath, 'spm12');
    
    addpath(spm_path);
    addpath(hmri_path);
    disp("path addded")

    run(path_to_batch)
    disp(matlabbatch)
    spm_jobman('run',matlabbatch)
end
