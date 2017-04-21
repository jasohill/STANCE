function [V_Tissue_Labels_Resliced,Y_Tissue_Labels_Resliced] = STANCE_reslice_tissue(fileName,scan,sumThreshold,gzipFlag,showFlag,This_sss)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% STANCE_reslice_tissue generates tissue fuzzy memberships from fileName in dataDir (optional) 
%   into functional space, that is specified according to the fMRI scan protocol,
%   resulting in a data structure Tissue_Labels_Resliced:
%
%        V = same format as SPM12 header
%        Y = MRI data array
%
%        scan.voxel.size: 3x1 list of voxel physical dimensions [mm]
%                        (default) = [3 3 3]
%        scan.voxel.matrix: 3x1 list of the number of voxels per spatial direction
%                        (default) = [64 64 NaN]
%                        NOTE: matrix(3) = NaN used indicate auto
%                              determine setting
%        scan.voxel.spacing: 3x1 list of inter-voxel spacing [mm]
%        scan.tiltAngle: tilt angle in YZ plane away from nose [degrees]
%                        (default) = 15
%
%        sumThreshold: (option 1) numerical threshold required upon a slice 
%                                 for inclusion in the resulting 3D volume.
%                      (option 2) a 2 element array where
%                                 sumThreshold(1) = the lower slice index
%                                 to use after the affine transformation
%                                 sumThreshold(1) = the upper slice index
%                                 to use after the affine transformation                             
%
%   NOTE: - if ANALYZE format, fileName should be *.hdr and *.img will be assumed to exist 
%         - can accomodate gziped file formats as well.
%         - optional mask for loading images can be provided.
%
%   Supported filetypes: NIFTI, ANALYZE and gzipped files
%
% Jason E. Hill
% STANCE_reslice_tissue.m      updated     13 JAN 2016
%------------------------------------------------------------------------
% 20 Anatomical Models of 20 Normal Brains
% "An new improved version of the realistic digital brain phantom" 
% - Berengere Aubert-Broche, Alan C. Evans, & Louis Collins
% NeuroImage 32 (2006) 138 - 145
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
% NOTE: effective T1's for background should be INF?
% * given as all 0. Values come from fitting to data.

% close all;

% global SPMversion 
% SPMversion = 'spm8';
% 
% %% add SPM to path (assumes toolbox is in the same directory as SPM)
% addpath(genpath(pwd));
% oldFolder = cd('../');
% cd(SPMversion);
% spmFolder = pwd;
% addpath(genpath(spmFolder));
% cd(oldFolder);

%% Turn off finite warning
warning('off', 'MATLAB:FINITE:obsoleteFunction');
%% Turn off nifti class warning when loading SPM mat files
warning('off', 'MATLAB:unknownElementsNowStruc');
warning('off', 'MATLAB:dispatcher:ShadowedMEXExtension');
warning('off', 'MATLAB:pfileOlderThanMfile');

%% Input argument handling and initialization

if nargin < 5
    showFlag = false;
end
if nargin < 4 
    gzipFlag = false;
end
if isempty(gzipFlag)
    gzipFlag = false;
end
if nargin < 3
    sumThreshold = 100;
end
if isempty(sumThreshold)
    sumThreshold = 100; 
end
if nargin < 2
    scan = [];
end

if ~exist('STANCE.mat','file')
    if ~exist('../STANCE.mat','file')
       load('../../STANCE.mat'); 
    else
       load('../STANCE.mat');
    end
else
    load('STANCE.mat');
end

if nargin < 6
    This_sss = Now_sss;
end

if nargin < 2 || isempty(scan)
    % assume defaults
    scan.voxel.size    = [3 3 3];
    scan.voxel.matrix  = [64 64 NaN]; 
    scan.voxel.spacing = [0 0 0.2*scan.voxel.size(3)]; % assume 20% Z-spacing
    scan.tiltAngle     = 15; % degrees 
end
new_dims = scan.voxel.matrix;
if isnan(new_dims(3))
    scan.sliceNbrFlag = false;
else
    sliceNbr = new_dims(3);
    scan.sliceNbrFlag = true;
end

%% Load the tissue labels
[V_Tissue_Labels,Y_Tissue_Labels] = STANCE_load_volume(fileName,[],gzipFlag);
if isempty(V_Tissue_Labels)
    error('Tissue label file failed to load.');
end

msg = 'Reslicing gray matter fuzzy membership labels...';
disp(msg);

% reslice gray matter fuzzy membership labels
Y_Gray_Matter = squeeze(Y_Tissue_Labels(:,:,:,3));
V_Gray_Matter = V_Tissue_Labels(3);
V_Gray_Matter.fname = 'gray_matter_anat.nii';
V_Gray_Matter.pinfo = [1 0 0]';
V_Gray_Matter.n = [];
V_Gray_Matter = spm_create_vol(V_Gray_Matter);
V_Gray_Matter = spm_write_vol(V_Gray_Matter,Y_Gray_Matter);

[V_Gray_Matter_Resliced,Y_Gray_Matter_Resliced] = STANCE_reslice_volume(V_Gray_Matter,scan,sumThreshold);

oldFilename = V_Gray_Matter_Resliced.fname;
V_Gray_Matter_Resliced.fname = [STANCE_genpath(This_sss,2),'/gray_matter_mask.nii'];
V_Gray_Matter_Resliced.pinfo = [1 0 0]';
V_Gray_Matter_Resliced.n = [];
%V_Gray_Matter_Resliced = spm_create_vol(V_Gray_Matter_Resliced);
spm_write_vol(V_Gray_Matter_Resliced,Y_Gray_Matter_Resliced);

sliceLimitLower = V_Gray_Matter_Resliced.sliceLimitLower;
sliceLimitUpper = V_Gray_Matter_Resliced.sliceLimitUpper;
sliceLimits = [sliceLimitLower, sliceLimitUpper];
showSlice = round(size(Y_Gray_Matter_Resliced,3)*0.6);

if showFlag
    figure; imshow(Y_Gray_Matter_Resliced(:,:,showSlice))
end

for k = 1:12
    V_Tissue_Labels_Resliced(k) = V_Gray_Matter_Resliced; %#ok<AGROW>
end
Y_Tissue_Labels_Resliced(:,:,:,3) = Y_Gray_Matter_Resliced;

delete(oldFilename)
delete('gray_matter_anat.nii')
clear Y_Gray_Matter;
clear V_Gray_Matter;
clear Y_Gray_Matter_Resliced;


fname = [STANCE_genpath(This_sss,2),'\background_mask.nii'];
msg = 'Reslicing background fuzzy membership labels...';
disp(msg);

% reslice background fuzzy membership labels
Y_Background       = squeeze(Y_Tissue_Labels(:,:,:,1));
V_Background       = V_Tissue_Labels(1);
V_Background.fname = 'background_anat.nii';
V_Background.pinfo = [1 0 0]';
V_Background.n     = [];
V_Background       = spm_create_vol(V_Background);
V_Background       = spm_write_vol(V_Background,Y_Background);

[V_Background_Resliced,Y_Background_Resliced] = STANCE_reslice_volume(V_Background,scan,sliceLimits);
oldFilename = V_Background_Resliced.fname;
V_Background_Resliced.fname = [STANCE_genpath(This_sss,2),'\background_mask.nii'];
V_Background_Resliced.pinfo = [1 0 0]';
V_Background_Resliced.n     = [];
V_Background_Resliced       = spm_create_vol(V_Background_Resliced);
spm_write_vol(V_Background_Resliced,Y_Background_Resliced);

Y_Tissue_Labels_Resliced(:,:,:,1) = Y_Background_Resliced;

if showFlag
figure; imshow(Y_Background_Resliced(:,:,showSlice))
end

delete(oldFilename);
delete('background_anat.nii');
clear Y_Background;
clear V_Background;
clear Y_Background_Resliced;
clear V_Background_Resliced;


msg = 'Reslicing CSF fuzzy membership labels...';
disp(msg);

% reslice CSF fuzzy membership labels
Y_CSF = squeeze(Y_Tissue_Labels(:,:,:,2));
V_CSF = V_Tissue_Labels(2);
V_CSF.fname = 'csf_anat.nii';
V_CSF.pinfo = [1 0 0]';
V_CSF.n = [];
V_CSF = spm_create_vol(V_CSF);
V_CSF = spm_write_vol(V_CSF,Y_CSF);

[V_CSF_Resliced,Y_CSF_Resliced] = STANCE_reslice_volume(V_CSF,scan,sliceLimits);
oldFilename = V_CSF_Resliced.fname;
V_CSF_Resliced.fname = [STANCE_genpath(This_sss,2),'\CSF_mask.nii'];
V_CSF_Resliced.pinfo = [1 0 0]';
V_CSF_Resliced.n = [];
V_CSF_Resliced = spm_create_vol(V_CSF_Resliced);
spm_write_vol(V_CSF_Resliced,Y_CSF_Resliced);

Y_Tissue_Labels_Resliced(:,:,:,2) = Y_CSF_Resliced;

if showFlag
    figure; imshow(Y_CSF_Resliced(:,:,showSlice))
end

delete(oldFilename);
delete('csf_anat.nii');
clear Y_CSF;
clear V_CSF;
clear Y_CSF_Resliced;
clear V_CSF_Resliced;


msg = 'Reslicing white matter fuzzy membership labels...';
disp(msg);

% reslice white matter fuzzy membership labels
Y_White_Matter = squeeze(Y_Tissue_Labels(:,:,:,4));
V_White_Matter = V_Tissue_Labels(4);
V_White_Matter.fname = 'white_matter_anat.nii';
V_White_Matter.pinfo = [1 0 0]';
V_White_Matter.n = [];
V_White_Matter = spm_create_vol(V_White_Matter);
V_White_Matter = spm_write_vol(V_White_Matter,Y_White_Matter);

[V_White_Matter_Resliced,Y_White_Matter_Resliced] = STANCE_reslice_volume(V_White_Matter,scan,sliceLimits);
oldFilename = V_White_Matter_Resliced.fname;

Y_Tissue_Labels_Resliced(:,:,:,4) = Y_White_Matter_Resliced;

if showFlag
    figure; imshow(Y_White_Matter_Resliced(:,:,showSlice))
end

delete(oldFilename);
delete('white_matter_anat.nii');
clear Y_White_Matter;
clear V_White_Matter;
clear Y_White_Matter_Resliced;
clear V_White_Matter_Resliced;


%% Fat tissue ...
msg = 'Reslicing fat fuzzy membership labels...';
disp(msg);

% reslice fat fuzzy membership labels
Y_Fat = squeeze(Y_Tissue_Labels(:,:,:,5));
V_Fat = V_Tissue_Labels(5);
V_Fat.fname = 'fat_anat.nii';
V_Fat.pinfo = [1 0 0]';
V_Fat.n = [];
V_Fat = spm_create_vol(V_Fat);
V_Fat = spm_write_vol(V_Fat,Y_Fat);

[V_Fat_Resliced,Y_Fat_Resliced] = STANCE_reslice_volume(V_Fat,scan,sliceLimits);
oldFilename = V_Fat_Resliced.fname;

Y_Tissue_Labels_Resliced(:,:,:,5) = Y_Fat_Resliced;

if showFlag
    figure; imshow(Y_Fat_Resliced(:,:,showSlice))
end

delete(oldFilename);
delete('fat_anat.nii');
clear Y_Fat;
clear V_Fat;
clear Y_Fat_Resliced;
clear V_Fat_Resliced;


%% Muscle tissue ...
msg = 'Reslicing muscle fuzzy membership labels...';
disp(msg);

% reslice muscle fuzzy membership labels
Y_Muscle = squeeze(Y_Tissue_Labels(:,:,:,6));
V_Muscle = V_Tissue_Labels(6);
V_Muscle.fname = 'muscle_anat.nii';
V_Muscle.pinfo = [1 0 0]';
V_Muscle.n = [];
V_Muscle = spm_create_vol(V_Muscle);
V_Muscle = spm_write_vol(V_Muscle,Y_Muscle);

[V_Muscle_Resliced,Y_Muscle_Resliced] = STANCE_reslice_volume(V_Muscle,scan,sliceLimits);
oldFilename = V_Muscle_Resliced.fname;

Y_Tissue_Labels_Resliced(:,:,:,6) = Y_Muscle_Resliced;

if showFlag
    figure; imshow(Y_Muscle_Resliced(:,:,showSlice))
    figure; imshow(Y_Muscle_Resliced(:,:,1))
end

delete(oldFilename);
delete('muscle_anat.nii');
clear Y_Muscle;
clear V_Muscle;
clear Y_Muscle_Resliced;
clear V_Muscle_Resliced;


%% Skin tissue ...
msg = 'Reslicing skin fuzzy membership labels...';
disp(msg);

% reslice skin fuzzy membership labels
Y_Skin = squeeze(Y_Tissue_Labels(:,:,:,7));
V_Skin = V_Tissue_Labels(7);
V_Skin.fname = 'skin_anat.nii';
V_Skin.pinfo = [1 0 0]';
V_Skin.n = [];
V_Skin = spm_create_vol(V_Skin);
V_Skin = spm_write_vol(V_Skin,Y_Skin);

[V_Skin_Resliced,Y_Skin_Resliced] = STANCE_reslice_volume(V_Skin,scan,sliceLimits);
oldFilename = V_Skin_Resliced.fname;

Y_Tissue_Labels_Resliced(:,:,:,7) = Y_Skin_Resliced;

if showFlag
    figure; imshow(Y_Skin_Resliced(:,:,showSlice))
end
 
Y_Skin_Resliced = imerode(Y_Skin_Resliced,ones(4,4,4));
Y_Skin_Resliced = imdilate(Y_Skin_Resliced,ones(4,4,4));

V_Skin_Resliced.fname = [STANCE_genpath(This_sss,2),'\eye_mask.nii'];
V_Skin_Resliced.pinfo = [1 0 0]';
V_Skin_Resliced.n = [];
V_Skin_Resliced = spm_create_vol(V_Skin_Resliced);
spm_write_vol(V_Skin_Resliced,Y_Skin_Resliced);

if showFlag
    figure; imshow(Y_Skin_Resliced(:,:,1))
end

delete(oldFilename);
delete('skin_anat.nii');
clear Y_Skin;
clear V_Skin;
clear Y_Skin_Resliced;
clear V_Skin_Resliced;


%% Skull tissue ...
msg = 'Reslicing skull fuzzy membership labels...';
disp(msg);

% reslice skull fuzzy membership labels
Y_Skull = squeeze(Y_Tissue_Labels(:,:,:,8));
V_Skull = V_Tissue_Labels(8);
V_Skull.fname = 'skull_anat.nii';
V_Skull.pinfo = [1 0 0]';
V_Skull.n = [];
V_Skull = spm_create_vol(V_Skull);
V_Skull = spm_write_vol(V_Skull,Y_Skull);

[V_Skull_Resliced,Y_Skull_Resliced] = STANCE_reslice_volume(V_Skull,scan,sliceLimits);
oldFilename = V_Skull_Resliced.fname;

Y_Tissue_Labels_Resliced(:,:,:,8) = Y_Skull_Resliced;

if showFlag
    figure; imshow(Y_Skull_Resliced(:,:,showSlice))
end

delete(oldFilename);
delete('skull_anat.nii');
clear Y_Skull;
clear V_Skull;
clear Y_Skull_Resliced;
clear V_Skull_Resliced;


% Blood vessels ...
msg = 'Reslicing blood vessels fuzzy membership labels...';
disp(msg);

% reslice vessels fuzzy membership labels
Y_Vessels = squeeze(Y_Tissue_Labels(:,:,:,9));
V_Vessels = V_Tissue_Labels(9);
V_Vessels.fname = 'vessels_anat.nii';
V_Vessels.pinfo = [1 0 0]';
V_Vessels.n = [];
V_Vessels = spm_create_vol(V_Vessels);
V_Vessels = spm_write_vol(V_Vessels,Y_Vessels);

[V_Vessels_Resliced,Y_Vessels_Resliced] = STANCE_reslice_volume(V_Vessels,scan,sliceLimits);
oldFilename = V_Vessels_Resliced.fname;
V_Vessels_Resliced.fname = [STANCE_genpath(This_sss,2),'\Vessels_mask.nii'];
V_Vessels_Resliced.pinfo = [1 0 0]';
V_Vessels_Resliced.n = [];
V_Vessels_Resliced = spm_create_vol(V_Vessels_Resliced);
spm_write_vol(V_Vessels_Resliced,Y_Vessels_Resliced);

Y_Tissue_Labels_Resliced(:,:,:,9) = Y_Vessels_Resliced;

if showFlag
    figure; imshow(Y_Vessels_Resliced(:,:,showSlice))
end

delete(oldFilename);
delete('vessels_anat.nii');
clear Y_Vessels;
clear V_Vessels;
clear Y_Vessels_Resliced;
clear V_Vessels_Resliced;


%% Connective tissue ...
msg = 'Reslicing connective tissue fuzzy membership labels...';
disp(msg);

% reslice connective (fat2) fuzzy membership labels
Y_Connective = squeeze(Y_Tissue_Labels(:,:,:,10));
V_Connective = V_Tissue_Labels(10);
V_Connective.fname = 'connective_anat.nii';
V_Connective.pinfo = [1 0 0]';
V_Connective.n = [];
V_Connective = spm_create_vol(V_Connective);
V_Connective = spm_write_vol(V_Connective,Y_Connective);

[V_Connective_Resliced,Y_Connective_Resliced] = STANCE_reslice_volume(V_Connective,scan,sliceLimits);
oldFilename = V_Connective_Resliced.fname;

Y_Tissue_Labels_Resliced(:,:,:,10) = Y_Connective_Resliced;

if showFlag
    figure; imshow(Y_Connective_Resliced(:,:,showSlice))
end

delete(oldFilename);
delete('connective_anat.nii');
clear Y_Connective;
clear V_Connective;
clear Y_Connective_Resliced;
clear V_Connective_Resliced;


%% Dura matter ...
msg = 'Reslicing dura matter fuzzy membership labels...';
disp(msg);

% reslice dura fuzzy membership labels
Y_Dura = squeeze(Y_Tissue_Labels(:,:,:,11));
V_Dura = V_Tissue_Labels(11);
V_Dura.fname = 'dura_anat.nii';
V_Dura.pinfo = [1 0 0]';
V_Dura.n = [];
V_Dura = spm_create_vol(V_Dura);
V_Dura = spm_write_vol(V_Dura,Y_Dura);

[V_Dura_Resliced,Y_Dura_Resliced] = STANCE_reslice_volume(V_Dura,scan,sliceLimits);
oldFilename = V_Dura_Resliced.fname;

Y_Tissue_Labels_Resliced(:,:,:,11) = Y_Dura_Resliced;

if showFlag
    figure; imshow(Y_Dura_Resliced(:,:,showSlice))
end

delete(oldFilename);
delete('dura_anat.nii');
clear Y_Dura;
clear V_Dura;
clear Y_Dura_Resliced;
clear V_Dura_Resliced;


%% Bone marrow ...
msg = 'Reslicing bone marrow fuzzy membership labels...';
disp(msg);

% reslice bone marrow fuzzy membership labels
Y_Marrow = squeeze(Y_Tissue_Labels(:,:,:,12));
V_Marrow = V_Tissue_Labels(12);
V_Marrow.fname = 'marrow_anat.nii';
V_Marrow.pinfo = [1 0 0]';
V_Marrow.n = [];
V_Marrow = spm_create_vol(V_Marrow);
V_Marrow = spm_write_vol(V_Marrow,Y_Marrow);

[V_Marrow_Resliced,Y_Marrow_Resliced] = STANCE_reslice_volume(V_Marrow,scan,sliceLimits);
oldFilename = V_Marrow_Resliced.fname;

Y_Tissue_Labels_Resliced(:,:,:,12) = Y_Marrow_Resliced;

if showFlag
    figure; imshow(Y_Marrow_Resliced(:,:,showSlice))
end

delete(oldFilename);
delete('marrow_anat.nii');
clear Y_Marrow;
clear V_Marrow;
clear Y_Marrow_Resliced;
clear V_Marrow_Resliced;

%% save resulting resliced tissue fuzzy membership data

Y_Tissue_Labels_Resliced(isnan(Y_Tissue_Labels_Resliced(:))) = 0;

fileExtension = fileName(end-2:end);
if strcmpi(fileExtension,'.gz')
    fileExtension = fileName(end-6:end-3);
    fileBase = fileName(1:end-7);
else
    fileBase = fileName(1:end-3);
end
fileBase = fileBase(end-14:end);

dataFileNameNIIout = [STANCE_genpath(This_sss,2),'\',fileBase,'_resliced',fileExtension]; 

V_Tissue_Labels_Resliced(1).fname = dataFileNameNIIout;
V_Tissue_Labels_Resliced_out      = V_Tissue_Labels_Resliced(1);
for k = 1:12
    V_Tissue_Labels_Resliced(k).fname = dataFileNameNIIout;
    V_Tissue_Labels_Resliced_out.n(1) = k;
    V_Tissue_Labels_Resliced(k) = spm_write_vol(V_Tissue_Labels_Resliced_out,squeeze(Y_Tissue_Labels_Resliced(:,:,:,k)));    
end

if gzipFlag
    dataFileNameNIIold = dataFileNameNIIout;
    dataFileNameNIIout = gzip(dataFileNameNIIout);
    delete(dataFileNameNIIold);
    V_Tissue_Labels_Resliced_out(:).fname = dataFileNameNIIout;
end

end
 
