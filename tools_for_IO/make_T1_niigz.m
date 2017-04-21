% Generates gzipped NIFTI of the T1-w images derived from McGill's BrainWeb data
%------------------------------------------------------------------------
% Author: Dr. Jason E. Hill, post-doctoral fellow
% Institution: CNG at TTU in a working association with TTNI
% Date: 28 July 2015
%------------------------------------------------------------------------
% 20 Anatomical Models of 20 Normal Brains
% NAME              relative PD     T1[ms]      T2[ms]      T2*[ms]
% Background*       0    5.0e-4          7.3e-4      0.3         0.001
% CSF               1    1.0             2569.0      329.0       58.0
% Grey Matter       2    0.86            833.0       83.0        69.0
% White Matter      3    0.77            500.0       70.0        61.0
% Fat               4    1.0             350.0       70.0        58.0
% Muscle            5    1.0             900.0       47.0        30.0
% Muscle/Skin       6    1.0             2569.0      329.0       58.0
% Skull*            7    0.26            2580.0      314.8       97.5
% Vessels*          8    1.0             800.0       180.0       58.0
% Connective (fat2) 9    0.77            500.0       70.0        61.0
% Dura Mater       10    0.77            2569.0      329.0       58.0
% Bone Marrow      11    1.0             500.0       70.0        61.0
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
       datafileNameMNC = ['subject0', num2str(subjects(s)), '_t1w_p4.mnc']; 
       datafileNameNII = ['subject0', num2str(subjects(s)), '_t1w_p0.nii'];  
   else
       datafileNameMNC = ['subject', num2str(subjects(s)), '_t1w_p4.mnc'];
       datafileNameNII = ['subject', num2str(subjects(s)), '_t1w_p0.nii'];
   end
   
   [T1_volRaw,T1_scaninfo] = loadminc(datafileNameMNC);
   figure; imshow(T1_volRaw(:,:,showSlice), []), drawnow;   
   
   T1_vol = T1_volRaw(21:237,39:219,:);

   VolInfo = spm_vol_minc(datafileNameMNC,1);
   VolInfo.dim(1) = 181;
   VolInfo.dim(2) = 217;       
   VolInfo.mat(1,4) = VolInfo.mat(1,4) + 38;
   VolInfo.mat(2,4) = VolInfo.mat(2,4) + 20;
   voxel_size = [1 1 1];
   origin = [-VolInfo.mat(1,4) -VolInfo.mat(2,4) -VolInfo.mat(3,4)]; 
   datatype = 16;
   
%    figure; imshow(T1_vol(:,:,showSlice), []), drawnow;   
%    nii = make_nii(T1_vol, voxel_size, origin, datatype); 
%    
%    save_nii(nii,datafileNameNII);
%    gzip(datafileNameNII);
%    delete(datafileNameNII);
   
   % denoise with bilateral filter
   T1min = min(min(min(T1_vol)));
   T1max = max(max(max(T1_vol)));
   T1_vol = (T1_vol-T1min)/(T1max-T1min);
      
   T1_vol = bfilter3(T1_vol,1,[3 0.06]);
    
   if T1min > 0 
       T1_vol = T1_vol*(T1max-T1min) + T1min;
   else
       T1_vol = T1_vol*T1max;
   end
    
   for z = 1:181
       T1_vol2(:,:,z) = flipud(T1_vol(:,:,z)');
   end
   
   figure; imshow(T1_vol2(:,:,showSlice), []), drawnow;
 
   nii = make_nii(T1_vol2, voxel_size, origin, datatype); 
   
   save_nii(nii,datafileNameNII);
   gzip(datafileNameNII);
   delete(datafileNameNII);
   
end

