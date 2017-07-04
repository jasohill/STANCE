function noiseMap = STANCE_make_noise_map(fn_tissues,dilation,factor,tissues,This_sss)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%
% Create a spatially varying noise amplitude map based on tissue priors.
%
% Author: Jason E. Hill, Ph. D.
% STANCE_make_noise_map.m      updated     28 FEB 2016
%------------------------------------------------------------------------
% 20 Anatomical Models of 20 Normal Brains from
% "An new improved version of the realistic digital brain phantom" 
% - Berengere Aubert-Broche, Alan C. Evans, & Louis Collins
%   NeuroImage 32 (2006) 138 - 145
%------------------------------------------------------------------------
% NAME              relative PD     T1[ms]      T2[ms]      T2*[ms]
% Background*       0    5.0e-4          7.3e-4      0.3         0.001
% CSF               1    1.0             2569.0      329.0       58.0
% Grey Matter       2    0.86            833.0       83.0        69.0
% White Matter      3    0.77            500.0       70.0        61.0
% Fat               4    1.0             350.0       70.0        58.0
% Muscle            5    1.0             900.0       47.0       *40.0
% Muscle/Skin       6    1.0            *1800.0     *140.0      *48.0
% Skull*            7    0.28            2580.0      314.8       50.0
% Vessels*          8    0.87            900.0       180.0       58.0
% Connective (fat2) 9    0.77            500.0       70.0        61.0
% Dura Mater       10    1.0             2569.0      329.0       58.0
% Bone Marrow      11    1.0            *2500.0      70.0        61.0

% handle variable input arguments
if nargin < 3
    factor = 2.0;
end
if factor < 1.0
    factor = 1.0 - factor;
end    
if nargin < 4
    if ~exist('STANCE.mat','file')
        if ~exist('..\STANCE.mat','file')
           load('..\..\STANCE.mat'); 
        else
           load('..\STANCE.mat');
        end
    else
        load('STANCE.mat');
    end
    if nargin < 5
       This_sss = Now_sss;
    end
    if isempty(fn_tissues)
        fileBase    = [STANCE_genpath(This_sss,2),'/'];
        CSFmask     = STANCE_load_data([fileBase,'CSF_mask.nii']);
        VesselsMask = STANCE_load_data([fileBase,'Vessels_mask.nii']);
        EyeMask     = STANCE_load_data([fileBase,'eye_mask.nii']);
    else
        % default: CSF and blood vessels + muscle/skin (surrogate for eyes)
        [~,Y_tissues] = STANCE_load_volume(fn_tissues);
        CSFmask     = Y_tissues(:,:,:,2);
        VesselsMask = Y_tissues(:,:,:,9);
        EyeMask     = Y_tissues(:,:,:,7);   
    end
    tissueMask = CSFmask + VesselsMask + 2*EyeMask;
else
    %e.g.: tissues = [2,9,7,7]; (CSF and blood vessels + eyes*2)
    [~,Y_tissues] = STANCE_load_volume(fn_tissues);  
    tissueMask = 0.*Y_tissues(:,:,:,tissues(1));
    for t = 1:length(tissues)
        tissueMask = tissueMask + Y_tissues(:,:,:,tissues(t));
    end
end
if dilation > 0
    noiseMap = (factor-1)*imdilate(tissueMask,ones(dilation,dilation,dilation)) + 1;
else
    noiseMap = (factor-1)*tissueMask + 1;    
end

