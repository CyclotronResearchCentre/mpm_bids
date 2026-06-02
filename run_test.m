paths = [%"/Volumes/ext4GB/Data/MPM/derivatives/BVSS/sub-C03/ses-20240812/mpm_cp/" 
    "/Volumes/ext4GB/Data/MPM/derivatives/BVSS/sub-C03/ses-20240812/mpm_ptx/"
    %"/Volumes/ext4GB/Data/MPM/derivatives/BVSS/sub-C04/ses-20240808/mpm_cp/"
    "/Volumes/ext4GB/Data/MPM/derivatives/BVSS/sub-C04/ses-20240808/mpm_ptx/"
    %"/Volumes/ext4GB/Data/MPM/derivatives/BVSS/sub-C10/ses-20240902/mpm_cp/"
    "/Volumes/ext4GB/Data/MPM/derivatives/BVSS/sub-C10/ses-20240902/mpm_ptx/"
    %"/Volumes/ext4GB/Data/MPM/derivatives/BVSS/sub-C11/ses-20240808/mpm_cp/"
    "/Volumes/ext4GB/Data/MPM/derivatives/BVSS/sub-C11/ses-20240808/mpm_ptx/"
    %"/Volumes/ext4GB/Data/MPM/derivatives/BVSS/sub-C14/ses-20240819/mpm_cp/"
    "/Volumes/ext4GB/Data/MPM/derivatives/BVSS/sub-C14/ses-20240819/mpm_ptx/"
    ];

for p = 1:1:14
    path = paths{p};
    list_batch = [fullfile(path,"spm_batch_helms_0p0.m") %fullfile(path,"spm_batch_helms_0p1.m")
        %fullfile(path,"spm_batch_helms_0p2.m") fullfile(path,"spm_batch_helms_0p3.m")
        %fullfile(path,"spm_batch_helms_0p4.m") fullfile(path,"spm_batch_helms_0p5.m")
        %fullfile(path,"spm_batch_helms_0p6.m") fullfile(path,"spm_batch_helms_0p7.m")
        %fullfile(path,"spm_batch_1p0.m") fullfile(path,"spm_batch_1p2.m") ...
        %fullfile(path,"spm_batch_1p4.m") fullfile(path,"spm_batch_1p6.m") ...
        %fullfile(path,"spm_batch_1p8.m") fullfile(path,"spm_batch_2p0.m") ...
        %fullfile(path,"spm_batch_2p2.m") fullfile(path,"spm_batch_2p4.m") ...
        ];
    disp(list_batch)
    for i = 1:1:length(list_batch)
        disp(list_batch{i})
        run(list_batch{i})
        disp(matlabbatch)
        spm_jobman('run',matlabbatch)
    end
end