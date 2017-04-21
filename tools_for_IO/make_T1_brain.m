% Generates gzipped NIFTI of the brain extracted from T1-w image derived from McGill's BrainWeb data
%------------------------------------------------------------------------
% Author: Dr. Jason E. Hill, post-doctoral fellow
% Institution: CNG at TTU in a working association with TTNI
% Date: 28 July 2015
%------------------------------------------------------------------------
% 20 Anatomical Models of 20 Normal Brains
% NAME              relative PD     T1[ms]      T2[ms]      T2*[ms]
%------------------------------------------------------------------.
% Background*       0    5.0e-4     7.3e-4      0.3         0.001  |
% CSF               1    1.0        2569.0      329.0       58.0   |
% Grey Matter       2    0.86       833.0       83.0        69.0   |
% White Matter      3    0.77       500.0       70.0        61.0   |
% Fat               4    1.0        350.0       70.0        58.0   |
% Muscle            5    1.0        900.0       47.0        30.0   |
% Muscle/Skin       6    1.0        2569.0      329.0       58.0   |
% Skull*            7    0.26       2580.0      314.8       97.5   | 
% Vessels*          8    1.0        800.0       180.0       58.0   |
% Connective (fat2) 9    0.77       500.0       70.0        61.0   |
% Dura Mater       10    1.0        2569.0      329.0       58.0   |
% Bone Marrow      11    1.0        500.0       70.0        61.0   |
%------------------------------------------------------------------'
% NOTE: effective T1's for background should be INF?
% * given as all 0.

clear all, close all;

addpath(genpath(pwd))

showSlice = 70; 

% %% LOADING RAW DATA
% %---------------------------------------
% % Load MINC raw data (source: BrainWeb)
% %---------------------------------------

subjects = [4,5,6,18,20,38,41,42,43,44,...
            45,46,47,48,49,50,51,52,53,54];

for s = 1:20
    
   if s < 4
       datafileNameMAT    = ['FuzzyBrain_s0', num2str(subjects(s)), '.mat'];
       datafileNameNIIin  = ['subject0', num2str(subjects(s)), '_t1w_p0.nii.gz'];
       datafileNameNIIout = ['subject0', num2str(subjects(s)), '_t1w_p0_brain.nii'];
   else
       datafileNameMAT    = ['FuzzyBrain_s', num2str(subjects(s)), '.mat'];
       datafileNameNIIin  = ['subject', num2str(subjects(s)), '_t1w_p0.nii.gz'];       
       datafileNameNIIout = ['subject', num2str(subjects(s)), '_t1w_p0_brain.nii'];       
   end
   
   FuzzyBrain = load(datafileNameMAT);
   nii = load_nii(datafileNameNIIin);
   
   [M,I] = max(FuzzyBrain.FuzzyBrain,[],4);
   
   for z = 1:181
       I2(:,:,z) = flipud(I(:,:,z)');
   end   
   
   % Remove skull and surrounding soft tissue: 
   %     brain = WM + GM (as compared with MNI152 brain)
   brainMask = (I2 == 3)|(I2 == 4);
   for z = 1:181
%      brainMask(:,:,z) = bwmorph(imopen(imfill(brainMask(:,:,z),'holes'),true(5)),'clean');
      brainMask(:,:,z) = bwmorph(imopen(brainMask(:,:,z),true(5)),'clean');
   end
   figure; imshow(brainMask(:,:,showSlice), []), drawnow;
   
   nii.img(~brainMask) = 0;
    
   figure; imshow(nii.img(:,:,showSlice), []), drawnow;
   
   save_nii(nii,datafileNameNIIout);
   gzip(datafileNameNIIout);
   delete(datafileNameNIIout);
   
end

