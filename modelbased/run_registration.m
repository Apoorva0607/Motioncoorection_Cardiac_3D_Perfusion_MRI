
function run_registration(fullcinemri1,rangey,rangex,referenceFrame, numSkip,numPreContrast_bld,numPreContrast_tiss,sliceNum)
num_slices=1;
% Preproceesing before deformable registration
% performs normalised cross correlation method for rigid shifts and then
% creates model images for deformable registration by fitting to 2 compartment model.

% Input: buffer_cinemri.mat (full FOV image)
% Output: cinemri_curv_fit (model images)

% set registration parameters here
numSkip=0;

noi=1;  % total number of MBR iterations


    
        
    try
        temp=strcat('Output1/',int2str(sliceNum),'/cinemri_curv_fit1.mat');
        load(temp)
     catch
        upsample_flag=0;upsample_factor=2;
        MBR_main(fullcinemri1,rangex,rangey,referenceFrame,numSkip,numPreContrast_bld,numPreContrast_tiss,noi,sliceNum,upsample_flag,upsample_factor);
    end    
        


