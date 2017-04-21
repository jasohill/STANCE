function [V_Parameter_Map,Y_Parameter_Map] = STANCE_make_parameter_map(filename,parameter,gzipFlag)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%STANCE_make_parameter_map generates a parameter map from the tissue fuzzy 
%   memberships provided by fileName in dataDir (optional) 
%
%        filename = filename of fuzzy data source
%        V = same format as SPM12 header
%        Y = MRI data array
%
%        parameter: 'T1'
%                   'T2'
%                   'T2*' or 'T2star'
%                   'PD'
%                   'delta'      - chemical shift 
%                   'PETatt'     - PET attenuation
%                   'PETact'     - PET activity
%                   'SPECTatt'   - SPECT attenuation
%                   'SPECTact'   - SPECT activity
%
%        baseLength: length of the base string of fileName  
%
%   NOTE: - if ANALYZE format, fileName should be *.hdr and *.img will be assumed to exist 
%         - can accomodate gziped file formats as well.
%         - optional mask for loading images can be provided.
%
%   Supported filetypes: NIFTI, ANALYZE and gzipped files
%
% Jason E. Hill
% STANCE_make_fun_tissue.m      updated     13 JAN 2016
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
% NOTE: effective T1's for background should be INF?
% * given as all 0. Values come from fitting to data.

if nargin < 3
    gzipFlag = false; %default 
end
if nargin < 2
    parameter = 'T2*'; %default 
end

%% Load the tissue fuzzy membership labels
%------------------------------------------
[V_Tissue_Labels,Y_Tissue_Labels] = STANCE_load_volume(filename);

%% Define tissue parameters

PETattenuation = [0.0, 0.09853, 0.09853, 0.0956, 0.151108, 0.09853, 0.09853, 0.087718, 0.087718, 0.098731, 0.098731, 0.0956]; % [cm^-1]
PETactivity    = [0.0, 22990.0, 8450.0,     0.0,      0.0,     0.0,     0.0,   8450.0,   8450.0,   8450.0,   8450.0,    0.0]; % [Bq/cm^3]

SPECTattenuation = [0.0, 0.1551, 0.1551, 0.1508, 0.3222, 0.1551, 0.1551, 0.1394, 0.1394, 0.1553, 0.1553, 0.1508]; % [cm^-1]
SPECTactivity    = [0.0, 11764.0,8194.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0]; % [Bq/cm^3]

% average tissue parameters as provided from BrainWeb for 1.5 T initialization:
T2_tissueAve0     = [0.3,     329.0,   83.0,  70.0,  70.0,  47.0,  140.0,  314.8, 180.0,  70.0,  329.0,   70.0]; % [ms]
delta_cs          = [0.0,     0.0,    0.0,   0.0,    -3.5,   0.0,    0.0,    0.0,   0.0,   0.0,    0.0,   -3.5]; % [ppm] at 1.5 T  

% First iteration correction of values as described in
% "Twenty New Digital Brain Phantoms for Creation of Validation Image Data Bases"
% - Berengere Aubert-Broche, Mark Griffin, G. Bruce Pike, Alan C. Evans, & D. Louis Collins
%   IEEE Trans. on Medical Imaging, Vol. 25, No. 11, NOV 2006
% method: estimate values from (T1raw  + T1bf + 4.0*T1_0)/6.0;
% where T1raw is the raw T1-w anatomical volume with 4% noise
%       T1bf is the bilateral filtered output of T1raw
%       T1_0 is the estimated
%       MSE = 0.0023 (best performing found)
PD_tissueAve1     = [0.3337,    1.0000,    0.8493,    0.7457,    0.9352,    0.9969,    1.0000,    0.1775,    0.7873,    0.7691,    1.0000,    1.0000];
T1_tissueAve1     = 1.0e+03.*[0.0000,    2.4214,    0.8500,    0.5249,    0.3916,    0.8972,    1.7282,    2.6350,    1.0490,    0.5030,    2.4364,    2.4671];
T2star_tissueAve1 = [0.0537,   61.9764,   64.5388,   53.0067,   46.7860,   40.2966,   63.4155,   45.7685,   44.9011,   59.0336,   69.7980,   65.9437];

%% Define fuzzy memberships
%---------------------------

background  = Y_Tissue_Labels(:,:,:,1);
CSF         = Y_Tissue_Labels(:,:,:,2);
grayMatter  = Y_Tissue_Labels(:,:,:,3);
whiteMatter = Y_Tissue_Labels(:,:,:,4);
fat         = Y_Tissue_Labels(:,:,:,5);
muscle      = Y_Tissue_Labels(:,:,:,6);
skin        = Y_Tissue_Labels(:,:,:,7);    
skull       = Y_Tissue_Labels(:,:,:,8);
vessels     = Y_Tissue_Labels(:,:,:,9);
connective  = Y_Tissue_Labels(:,:,:,10);
dura        = Y_Tissue_Labels(:,:,:,11);
marrow      = Y_Tissue_Labels(:,:,:,12);

% showSlice = 20;
% figure,imshow(background(:,:,showSlice),[]);
% figure,imshow(CSF(:,:,showSlice),[]);
% figure,imshow(grayMatter(:,:,showSlice),[]);
% figure,imshow(whiteMatter(:,:,showSlice),[]);

%% Generate maps from fuzzy labels
%------------------------------------------
%                   'T2'
%                   'T2*'
%                   'PD'
%                   'delta'      - chemical shift 
%                   'PETatt'     - PET attenuation
%                   'PETact'     - PET activity
%                   'SPECTatt'   - SPECT attenuation
%                   'SPECTact'
V_Parameter_Map = V_Tissue_Labels(1);

if strcmp(parameter,'T1')
   Y_Parameter_Map = T1_tissueAve1(1)*background + T1_tissueAve1(2)*CSF + T1_tissueAve1(3)*grayMatter ...
       + T1_tissueAve1(4)*whiteMatter + T1_tissueAve1(5)*fat + T1_tissueAve1(6)*muscle + T1_tissueAve1(7)*skin ...
       + T1_tissueAve1(8)*skull + T1_tissueAve1(9)*vessels + T1_tissueAve1(10)*connective + T1_tissueAve1(11)*dura + T1_tissueAve1(12)*marrow;   
   Pmax = max(T1_tissueAve1);
   Pmin = 0.0;
end
if strcmp(parameter,'T2')
   Y_Parameter_Map = T2_tissueAve0(1)*background + T2_tissueAve0(2)*CSF + T2_tissueAve0(3)*grayMatter ...
       + T2_tissueAve0(4)*whiteMatter + T2_tissueAve0(5)*fat + T2_tissueAve0(6)*muscle + T2_tissueAve0(7)*skin ...
       + T2_tissueAve0(8)*skull + T2_tissueAve0(9)*vessels + T2_tissueAve0(10)*connective + T2_tissueAve0(11)*dura + T2_tissueAve0(12)*marrow;  
   Pmax = max(T2_tissueAve0);
   Pmin = 0.0;
end
if strcmp(parameter,'T2*') || strcmp(parameter,'T2star') 
   Y_Parameter_Map = T2star_tissueAve1(1)*background + T2star_tissueAve1(2)*CSF + T2star_tissueAve1(3)*grayMatter ...
       + T2star_tissueAve1(4)*whiteMatter + T2star_tissueAve1(5)*fat + T2star_tissueAve1(6)*muscle + T2star_tissueAve1(7)*skin ...
       + T2star_tissueAve1(8)*skull + T2star_tissueAve1(9)*vessels + T2star_tissueAve1(10)*connective + T2star_tissueAve1(11)*dura + T2star_tissueAve1(12)*marrow;  
   Pmax = max(T2star_tissueAve1);
   Pmin = 0.0;
end
if strcmp(parameter,'PD')   
   Y_Parameter_Map = PD_tissueAve1(1)*background + PD_tissueAve1(2)*CSF + PD_tissueAve1(3)*grayMatter ...
       + PD_tissueAve1(4)*whiteMatter + PD_tissueAve1(5)*fat + PD_tissueAve1(6)*muscle + PD_tissueAve1(7)*skin ...
       + PD_tissueAve1(8)*skull + PD_tissueAve1(9)*vessels + PD_tissueAve1(10)*connective + PD_tissueAve1(11)*dura + PD_tissueAve1(12)*marrow;
   Pmax = max(PD_tissueAve1);
   Pmin = 0.0;
end
if strcmp(parameter,'delta')
   Y_Parameter_Map = delta_cs(1)*background + delta_cs(2)*CSF + delta_cs(3)*grayMatter ...
       + delta_cs(4)*whiteMatter + delta_cs(5)*fat + delta_cs(6)*muscle + delta_cs(7)*skin ...
       + delta_cs(8)*skull + delta_cs(9)*vessels + delta_cs(10)*connective + delta_cs(11)*dura + delta_cs(12)*marrow;    
   Pmax = max(delta_cs);
   Pmin = min(delta_cs);
end
if strcmp(parameter,'PETatt')
   Y_Parameter_Map = PETattenuation(1)*background + PETattenuation(2)*CSF + PETattenuation(3)*grayMatter ...
       + PETattenuation(4)*whiteMatter + PETattenuation(5)*fat + PETattenuation(6)*muscle + PETattenuation(7)*skin ...
       + PETattenuation(8)*skull + PETattenuation(9)*vessels + PETattenuation(10)*connective + PETattenuation(11)*dura + PETattenuation(12)*marrow;    
   Pmax = max(PETattenuation);
   Pmin = 0.0;
end
if strcmp(parameter,'PETact')
   Y_Parameter_Map = PETactivity(1)*background + PETactivity(2)*CSF + PETactivityn(3)*grayMatter ...
       + PETactivity(4)*whiteMatter + PETactivity(5)*fat + PETactivity(6)*muscle + PETactivity(7)*skin ...
       + PETactivity(8)*skull + PETactivity(9)*vessels + PETactivity(10)*connective + PETactivity(11)*dura + PETactivity(12)*marrow;    
   Pmax = max(PETactivity);
   Pmin = 0.0;
end
if strcmp(parameter,'SPECTatt')
   Y_Parameter_Map = SPECTattenuation(1)*background + SPECTattenuation(2)*CSF + SPECTattenuation(3)*grayMatter ...
       + SPECTattenuation(4)*whiteMatter + SPECTattenuation(5)*fat + SPECTattenuation(6)*muscle + SPECTattenuation(7)*skin ...
       + SPECTattenuation(8)*skull + SPECTattenuation(9)*vessels + SPECTattenuation(10)*connective + SPECTattenuation(11)*dura + SPECTattenuation(12)*marrow;    
   Pmax = max(SPECTattenuation);
   Pmin = 0.0;
end
if strcmp(parameter,'SPECTact')
   Y_Parameter_Map = SPECTactivity(1)*background + SPECTactivity(2)*CSF + SPECTactivityn(3)*grayMatter ...
       + SPECTactivity(4)*whiteMatter + SPECTactivity(5)*fat + SPECTactivity(6)*muscle + SPECTactivity(7)*skin ...
       + SPECTactivity(8)*skull + SPECTactivity(9)*vessels + SPECTactivity(10)*connective + SPECTactivity(11)*dura + SPECTactivity(12)*marrow;    
   Pmax = max(SPECTactivity);
   Pmin = 0.0;
end

%% write output file

fileExtension = filename(end-2:end);
if strcmpi(fileExtension,'.gz')
    dataFileNameNIIin = gunzip(filename);
    fileBase = filename(1:end-7);
else
    dataFileNameNIIin = filename;
    fileBase = filename(1:end-4);  
end

if strcmp(parameter,'T2*')
    parameter = 'T2star';
end

fileBase = [fileBase,'_',parameter,'_map'];
dataFileNameNIIout = [fileBase,'.nii'];
V_Parameter_Map.fname = dataFileNameNIIout;

dataDir = [];
if isempty(dataDir)
    V_Parameter_Map = spm_write_vol(V_Parameter_Map,squeeze(Y_Parameter_Map));    
else
    filePath = fullfile(dataDir,filename);
    oldPath = cd(filePath);
    V_Parameter_Map = spm_write_vol(V_Parameter_Map,squeeze(Y_Parameter_Map));    
    cd(oldPath);
end

if gzipFlag
    dataFileNameNIIold = dataFileNameNIIout;
    dataFileNameNIIout = gzip(dataFileNameNIIout);
    delete(dataFileNameNIIold);
    V_Tissue_Labels_Resliced_out(:).fname = dataFileNameNIIout;
end
   
end

