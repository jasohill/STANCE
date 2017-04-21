% Generates gzipped NIFTI derived from McGill's BrainWeb data
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
% Muscle            5    1.0             900.0       47.0       *40.0
% Muscle/Skin       6    1.0            *1800.0     *140.0      *48.0
% Skull*            7    0.28            2580.0      314.8       50.0
% Vessels*          8    0.87            900.0       180.0       58.0
% Connective (fat2) 9    0.77            500.0       70.0        61.0
% Dura Mater       10    1.0             2569.0      329.0       58.0
% Bone Marrow      11    1.0            *2500.0      70.0        61.0
% NOTE: effective T1's for background should be INF?
% * given as all 0.

clear all, close all;

addpath(genpath(pwd))

showSlice = 70; 

write2file = false;  

% Average tissue parameters as provided from BrainWeb for 1.5 T
% initialization
PD_tissueAve0     = [5.0e-4,    1.0,   0.86,  0.77,   1.0,   1.0,    1.0,   0.18,   0.87,  0.77,    1.0,   1.0]; % [rel. to H20 = 1.0]
T1_tissueAve0     = [7.3e-4, 2569.0,  833.0, 500.0, 350.0, 900.0, 1800.0, 2580.0, 900.0, 500.0, 2569.0, 2500.0]; % [ms]
T2_tissueAve0     = [0.3,     329.0,   83.0,  70.0,  70.0,  47.0,  140.0,  314.8, 180.0,  70.0,  329.0,   70.0]; % [ms]
T2star_tissueAve0 = [0.001,    58.0,   69.0,  61.0,  58.0,  40.0,   48.0,   50.0,  58.0,  61.0,   58.0,   61.0]; % [ms]
% NOTE: background and skull values are extrapolated from data. 
% MSE = 0.0031

% method = (raw + bf + 4.0*0)/6.0;
PD_tissueAve1     = [0.3337,    1.0000,    0.8493,    0.7457,    0.9352,    0.9969,    1.0000,    0.1775,    0.7873,    0.7691,    1.0000,    1.0000];
T1_tissueAve1     = 1.0e+03.*[0.0000,    2.4214,    0.8500,    0.5249,    0.3916,    0.8972,    1.7282,    2.6350,    1.0490,    0.5030,    2.4364,    2.4671];
T2star_tissueAve1 = [0.0537,   61.9764,   64.5388,   53.0067,   46.7860,   40.2966,   63.4155,   45.7685,   44.9011,   59.0336,   69.7980,   65.9437];
% MSE = 0.0023

% second run results from 1 updated with 0
PD_tissueAve20     = [0.3337,    0.9999,    0.8566,    0.7602,    0.9687,    0.9965,    0.9980,    0.1796,    0.8265,    0.7718,    0.9996,    0.9996];
T1_tissueAve20     = 1.0e+03.*[0.0000,    2.4288,    0.8424,    0.5119,    0.3707,    0.8966,    1.7554,    2.6090,    0.9820,    0.5015,    2.4571,    2.4826];
T2star_tissueAve20 = [0.3353,   81.5041,   67.9182,   57.9628,   52.1308,   40.3271,   59.5641,   46.0669   50.4677,   59.3153,   79.6380,   72.3186];
% MSE = 0.0025

% second run results from 1 updated with 1
PD_tissueAve21     = [0.5558,    0.9999,    0.8495,    0.7440,    0.9255,    0.9945,    0.9980,    0.1780,    0.7714,    0.7712,    0.9996,    0.9996];
T1_tissueAve21     = 1.0e+03.*[0,    2.3304,    0.8538,    0.5285,    0.3984,    0.8947,    1.7075,    2.6457,    1.0813,    0.5035,    2.3687,    2.4606];
T2star_tissueAve21 = [0.3704,   84.1551,   64.9440,   52.6339,   44.6548,   40.5248,   69.8411,   43.2459,   41.7351,   58.0044,   87.5033,   75.6144];
% MSE = 0.0025

MSE0 = zeros(1,20);
MSE1 = zeros(1,20);
MSE20 = zeros(1,20);
MSE21 = zeros(1,20);

% %% LOADING RAW DATA
% %---------------------------------------
% % Load MINC raw data (source: BrainWeb)
% %---------------------------------------

subjects = [4,5,6,18,20,38,41,42,43,44,...
            45,46,47,48,49,50,51,52,53,54];      
        
for s = 1:20   
%s = 1

   %% LOADING DATA
   %---------------------------------------
   % Load MINC raw data (source: BrainWeb)
   %--------------------------------------- 
   if s < 4
       datafileNameMAT = ['FuzzyBrain_s0', num2str(subjects(s)), '.mat'];
       datafileNameMNC = ['subject0', num2str(subjects(s)), '_t1w_p4.mnc']; 
       datafileNameNII = ['subject0', num2str(subjects(s)), '_t1w_p0.nii.gz'];
       datafileNameMATout = ['ParameterMaps_s0', num2str(subjects(s)), '.mat'];
   else
       datafileNameMAT = ['FuzzyBrain_s', num2str(subjects(s)), '.mat'];
       datafileNameMNC = ['subject', num2str(subjects(s)), '_t1w_p4.mnc'];
       datafileNameNII = ['subject', num2str(subjects(s)), '_t1w_p0.nii.gz'];
       datafileNameMATout = ['ParameterMaps_s', num2str(subjects(s)), '.mat'];       
   end

   [T1_volRaw,T1scaninfo] = loadminc(datafileNameMNC);
   T1raw = T1_volRaw(21:237,39:219,:); % T1-w raw data with 4% noise 
   nii = load_nii(datafileNameNII);
   T1bf = nii.img;                     % T1-w data bilater filtered
   
   % Load fuzzy memberships
   %-------------------------------   
   
   FuzzyBrain = load(datafileNameMAT);

   % define fuzzy memberships
   background  = double(zeros(217,181,181));
   CSF         = background;
   grayMatter  = background;
   whiteMatter = background;
   fat         = background;
   muscle      = background;
   muskin      = background;
   skull       = background;
   vessels     = background;
   fat2        = background;
   dura        = background;   
   marrow      = background; 
   
   for k = 1:181
       background(:,:,k)  = FuzzyBrain.FuzzyBrain(:,:,k,1);
       CSF(:,:,k)         = FuzzyBrain.FuzzyBrain(:,:,k,2);
       grayMatter(:,:,k)  = FuzzyBrain.FuzzyBrain(:,:,k,3);
       whiteMatter(:,:,k) = FuzzyBrain.FuzzyBrain(:,:,k,4);
       fat(:,:,k)         = FuzzyBrain.FuzzyBrain(:,:,k,5);
       muscle(:,:,k)      = FuzzyBrain.FuzzyBrain(:,:,k,6);
       muskin(:,:,k)      = FuzzyBrain.FuzzyBrain(:,:,k,7);    
       skull(:,:,k)       = FuzzyBrain.FuzzyBrain(:,:,k,8);
       vessels(:,:,k)     = FuzzyBrain.FuzzyBrain(:,:,k,9);
       fat2(:,:,k)        = FuzzyBrain.FuzzyBrain(:,:,k,10);
       dura(:,:,k)        = FuzzyBrain.FuzzyBrain(:,:,k,11);
       marrow(:,:,k)      = FuzzyBrain.FuzzyBrain(:,:,k,12);    
   end   

   % define fuzzy 100% masks
   background_1  = (background  > 0.985);
   CSF_1         = (CSF         > 0.985);
   grayMatter_1  = (grayMatter  > 0.985);
   whiteMatter_1 = (whiteMatter > 0.985);
   fat_1         = (fat         > 0.985);
   muscle_1      = (muscle      > 0.985);
   muskin_1      = (muskin      > 0.985);    
   skull_1       = (skull       > 0.985); 
   vessels_1     = (vessels     > 0.985); 
   fat2_1        = (fat2        > 0.985);    
   dura_1        = (dura        > 0.985); 
   marrow_1      = (marrow      > 0.985);    
    
   %% TISSUE & DATA STATISTICS
   %---------------------------------------
   
   % get tissue statistics of raw data
%   T1_tissueAve
   T1wRawStd = [ std(T1raw(background_1))  , std(T1raw(CSF_1))  , std(T1raw(grayMatter_1))  , std(T1raw(whiteMatter_1))  , std(T1raw(fat_1))  , std(T1raw(muscle_1))  , std(T1raw(muskin_1))  , std(T1raw(skull_1))  , std(T1raw(vessels_1))  , std(T1raw(fat2_1))  , std(T1raw(dura_1))  , std(T1raw(dura_1))  ]
   T1wRawMin = [ min(T1raw(background_1))  , min(T1raw(CSF_1))  , min(T1raw(grayMatter_1))  , min(T1raw(whiteMatter_1))  , min(T1raw(fat_1))  , min(T1raw(muscle_1))  , min(T1raw(muskin_1))  , min(T1raw(skull_1))  , min(T1raw(vessels_1))  , min(T1raw(fat2_1))  , min(T1raw(dura_1))  , min(T1raw(dura_1))  ]; 
   T1wRawMed = [median(T1raw(background_1)),median(T1raw(CSF_1)),median(T1raw(grayMatter_1)),median(T1raw(whiteMatter_1)),median(T1raw(fat_1)),median(T1raw(muscle_1)),median(T1raw(muskin_1)),median(T1raw(skull_1)),median(T1raw(vessels_1)),median(T1raw(fat2_1)),median(T1raw(dura_1)),median(T1raw(dura_1))] 
   T1wRawMax = [ max(T1raw(background_1))  , max(T1raw(CSF_1))  , max(T1raw(grayMatter_1))  , max(T1raw(whiteMatter_1))  , max(T1raw(fat_1))  , max(T1raw(muscle_1))  , max(T1raw(muskin_1))  , max(T1raw(skull_1))  , max(T1raw(vessels_1))  , std(T1raw(fat2_1))  , std(T1raw(dura_1))  , std(T1raw(dura_1))  ]; 

   T1wBfStd = [ std(T1bf(background_1))  , std(T1bf(CSF_1))  , std(T1bf(grayMatter_1))  , std(T1bf(whiteMatter_1))  , std(T1bf(fat_1))  , std(T1bf(muscle_1))  , std(T1bf(muskin_1))  , std(T1bf(skull_1))  , std(T1bf(vessels_1))  , std(T1bf(fat2_1))  , std(T1bf(dura_1))  , std(T1bf(dura_1))  ] 
   T1wBfMin = [ min(T1bf(background_1))  , min(T1bf(CSF_1))  , min(T1bf(grayMatter_1))  , min(T1bf(whiteMatter_1))  , min(T1bf(fat_1))  , min(T1bf(muscle_1))  , min(T1bf(muskin_1))  , min(T1bf(skull_1))  , min(T1bf(vessels_1))  , min(T1bf(fat2_1))  , min(T1bf(dura_1))  , min(T1bf(dura_1))  ]; 
   T1wBfMed = [median(T1bf(background_1)),median(T1bf(CSF_1)),median(T1bf(grayMatter_1)),median(T1bf(whiteMatter_1)),median(T1bf(fat_1)),median(T1bf(muscle_1)),median(T1bf(muskin_1)),median(T1bf(skull_1)),median(T1bf(vessels_1)),median(T1bf(fat2_1)),median(T1bf(dura_1)),median(T1bf(dura_1))] 
   T1wBfMax = [ max(T1bf(background_1))  , max(T1bf(CSF_1))  , max(T1bf(grayMatter_1))  , max(T1bf(whiteMatter_1))  , max(T1bf(fat_1))  , max(T1bf(muscle_1))  , max(T1bf(muskin_1))  , max(T1bf(skull_1))  , max(T1bf(vessels_1))  , min(T1bf(fat2_1))  , min(T1bf(dura_1))  , min(T1bf(dura_1))  ];    
   
   showSlice = 70;   
   Iraw = T1raw(:,:,showSlice)/max(T1wRawMax);
   Ibf  = double(T1bf(:,:,showSlice)'/max(T1wRawMax));
   figure, imshow(Iraw);
   figure, imshow(Ibf);   

   %% Generate initial maps from fuzzy labels
   %------------------------------------------
   PD_map0 = PD_tissueAve0(1)*background + PD_tissueAve0(2)*CSF + PD_tissueAve0(3)*grayMatter ...
       + PD_tissueAve0(4)*whiteMatter + PD_tissueAve0(5)*fat + PD_tissueAve0(6)*muscle + PD_tissueAve0(7)*muskin ...
       + PD_tissueAve0(8)*skull + PD_tissueAve0(9)*vessels + PD_tissueAve0(10)*fat2 + PD_tissueAve0(11)*dura + PD_tissueAve0(12)*marrow;
   T1_map0 = T1_tissueAve0(1)*background + T1_tissueAve0(2)*CSF + T1_tissueAve0(3)*grayMatter ...
       + T1_tissueAve0(4)*whiteMatter + T1_tissueAve0(5)*fat + T1_tissueAve0(6)*muscle + T1_tissueAve0(7)*muskin ...
       + T1_tissueAve0(8)*skull + T1_tissueAve0(9)*vessels + T1_tissueAve0(10)*fat2 + T1_tissueAve0(11)*dura + T1_tissueAve0(12)*marrow;   
   T2_map0 = T2_tissueAve0(1)*background + T2_tissueAve0(2)*CSF + T2_tissueAve0(3)*grayMatter ...
       + T2_tissueAve0(4)*whiteMatter + T2_tissueAve0(5)*fat + T2_tissueAve0(6)*muscle + T2_tissueAve0(7)*muskin ...
       + T2_tissueAve0(8)*skull + T2_tissueAve0(9)*vessels + T2_tissueAve0(10)*fat2 + T2_tissueAve0(11)*dura + T2_tissueAve0(12)*marrow;  
   T2star_map0 = T2star_tissueAve0(1)*background + T2star_tissueAve0(2)*CSF + T2star_tissueAve0(3)*grayMatter ...
       + T2star_tissueAve0(4)*whiteMatter + T2star_tissueAve0(5)*fat + T2star_tissueAve0(6)*muscle + T2star_tissueAve0(7)*muskin ...
       + T2star_tissueAve0(8)*skull + T2star_tissueAve0(9)*vessels + T2star_tissueAve0(10)*fat2 + T2star_tissueAve0(11)*dura + T2star_tissueAve0(12)*marrow;  
   
   PD_map1 = PD_tissueAve1(1)*background + PD_tissueAve1(2)*CSF + PD_tissueAve1(3)*grayMatter ...
       + PD_tissueAve1(4)*whiteMatter + PD_tissueAve1(5)*fat + PD_tissueAve1(6)*muscle + PD_tissueAve1(7)*muskin ...
       + PD_tissueAve1(8)*skull + PD_tissueAve1(9)*vessels + PD_tissueAve1(10)*fat2 + PD_tissueAve1(11)*dura + PD_tissueAve1(12)*marrow;
   T1_map1 = T1_tissueAve1(1)*background + T1_tissueAve1(2)*CSF + T1_tissueAve1(3)*grayMatter ...
       + T1_tissueAve1(4)*whiteMatter + T1_tissueAve1(5)*fat + T1_tissueAve1(6)*muscle + T1_tissueAve1(7)*muskin ...
       + T1_tissueAve1(8)*skull + T1_tissueAve1(9)*vessels + T1_tissueAve1(10)*fat2 + T1_tissueAve1(11)*dura + T1_tissueAve1(12)*marrow;   
   T2star_map1 = T2star_tissueAve1(1)*background + T2star_tissueAve1(2)*CSF + T2star_tissueAve1(3)*grayMatter ...
       + T2star_tissueAve1(4)*whiteMatter + T2star_tissueAve1(5)*fat + T2star_tissueAve1(6)*muscle + T2star_tissueAve1(7)*muskin ...
       + T2star_tissueAve1(8)*skull + T2star_tissueAve1(9)*vessels + T2star_tissueAve1(10)*fat2 + T2star_tissueAve1(11)*dura + T2star_tissueAve1(12)*marrow;  

   PD_map20 = PD_tissueAve20(1)*background + PD_tissueAve20(2)*CSF + PD_tissueAve20(3)*grayMatter ...
       + PD_tissueAve20(4)*whiteMatter + PD_tissueAve20(5)*fat + PD_tissueAve20(6)*muscle + PD_tissueAve20(7)*muskin ...
       + PD_tissueAve20(8)*skull + PD_tissueAve20(9)*vessels + PD_tissueAve20(10)*fat2 + PD_tissueAve20(11)*dura + PD_tissueAve20(12)*marrow;
   T1_map20 = T1_tissueAve20(1)*background + T1_tissueAve20(2)*CSF + T1_tissueAve20(3)*grayMatter ...
       + T1_tissueAve20(4)*whiteMatter + T1_tissueAve20(5)*fat + T1_tissueAve20(6)*muscle + T1_tissueAve20(7)*muskin ...
       + T1_tissueAve20(8)*skull + T1_tissueAve20(9)*vessels + T1_tissueAve20(10)*fat2 + T1_tissueAve20(11)*dura + T1_tissueAve20(12)*marrow;   
   T2star_map20 = T2star_tissueAve20(1)*background + T2star_tissueAve20(2)*CSF + T2star_tissueAve20(3)*grayMatter ...
       + T2star_tissueAve20(4)*whiteMatter + T2star_tissueAve20(5)*fat + T2star_tissueAve20(6)*muscle + T2star_tissueAve20(7)*muskin ...
       + T2star_tissueAve20(8)*skull + T2star_tissueAve20(9)*vessels + T2star_tissueAve20(10)*fat2 + T2star_tissueAve20(11)*dura + T2star_tissueAve20(12)*marrow;       
   
   PD_map21 = PD_tissueAve21(1)*background + PD_tissueAve21(2)*CSF + PD_tissueAve21(3)*grayMatter ...
       + PD_tissueAve21(4)*whiteMatter + PD_tissueAve21(5)*fat + PD_tissueAve21(6)*muscle + PD_tissueAve21(7)*muskin ...
       + PD_tissueAve21(8)*skull + PD_tissueAve21(9)*vessels + PD_tissueAve21(10)*fat2 + PD_tissueAve21(11)*dura + PD_tissueAve21(12)*marrow;
   T1_map21 = T1_tissueAve21(1)*background + T1_tissueAve21(2)*CSF + T1_tissueAve21(3)*grayMatter ...
       + T1_tissueAve21(4)*whiteMatter + T1_tissueAve21(5)*fat + T1_tissueAve21(6)*muscle + T1_tissueAve21(7)*muskin ...
       + T1_tissueAve21(8)*skull + T1_tissueAve21(9)*vessels + T1_tissueAve21(10)*fat2 + T1_tissueAve21(11)*dura + T1_tissueAve21(12)*marrow;   
   T2star_map21 = T2star_tissueAve21(1)*background + T2star_tissueAve21(2)*CSF + T2star_tissueAve21(3)*grayMatter ...
       + T2star_tissueAve21(4)*whiteMatter + T2star_tissueAve21(5)*fat + T2star_tissueAve21(6)*muscle + T2star_tissueAve21(7)*muskin ...
       + T2star_tissueAve21(8)*skull + T2star_tissueAve21(9)*vessels + T2star_tissueAve21(10)*fat2 + T2star_tissueAve21(11)*dura + T2star_tissueAve21(12)*marrow;    
   
   %% Generate synthetic data from initial maps
   %---------------------------------------------
   Mz         = 1.0;  % magnetization factor in z direction
   TRge       = 22;        % [ms]
   TEge       = 9.2;        % [ms]
   alpha      = 30*pi/180; % [degrees]
   T1synthetic = sin(alpha)*Mz*PD_map0.*((1-exp(-TRge/T1_map0))./(1-cos(alpha)*exp(-TRge/T1_map0))).*exp(-TEge/T2star_map0);
 
   % determine the magnetization multiplier Mz
   Nbr_CSF = sum(sum(sum(CSF_1)));
   Nbr_GM  = sum(sum(sum(grayMatter_1)));
   Nbr_WM  = sum(sum(sum(whiteMatter_1)));
   
   Mz = (T1wRawMed(2)/median(T1synthetic(CSF_1))*Nbr_CSF+...
   + T1wRawMed(3)/median(T1synthetic(grayMatter_1))*Nbr_GM + ...
   + T1wRawMed(4)/median(T1synthetic(whiteMatter_1))*Nbr_WM)/(Nbr_CSF+Nbr_GM+Nbr_WM)

   T1synthetic = sin(alpha)*Mz*PD_map0.*((1-exp(-TRge/T1_map0))./(1-cos(alpha)*exp(-TRge/T1_map0))).*exp(-TEge/T2star_map0);
   T1synthetic1 = sin(alpha)*Mz*PD_map1.*((1-exp(-TRge/T1_map1))./(1-cos(alpha)*exp(-TRge/T1_map1))).*exp(-TEge/T2star_map1);
   T1synthetic20 = sin(alpha)*Mz*PD_map20.*((1-exp(-TRge/T1_map20))./(1-cos(alpha)*exp(-TRge/T1_map20))).*exp(-TEge/T2star_map20); 
   T1synthetic21 = sin(alpha)*Mz*PD_map21.*((1-exp(-TRge/T1_map21))./(1-cos(alpha)*exp(-TRge/T1_map21))).*exp(-TEge/T2star_map21); 
   Isyn0 = T1synthetic(:,:,showSlice)/max(T1wRawMax);
   figure, imshow(Isyn0);
   MSE0(s) = 0.5*(immse(Iraw,Isyn0)+immse(Ibf,Isyn0));
   Isyn1 = T1synthetic1(:,:,showSlice)/max(T1wRawMax);
   MSE1(s) = 0.5*(immse(Iraw,Isyn1)+immse(Ibf,Isyn1));
   Isyn20 = T1synthetic20(:,:,showSlice)/max(T1wRawMax);
   MSE20(s) = 0.5*(immse(Iraw,Isyn20)+immse(Ibf,Isyn20));   
   Isyn21 = T1synthetic21(:,:,showSlice)/max(T1wRawMax);
   MSE21(s) = 0.5*(immse(Iraw,Isyn21)+immse(Ibf,Isyn21));
   
   % validate
   CSFnormRaw = T1wRawMed(2)/median(T1synthetic(CSF_1))
   gmNormRaw  = T1wRawMed(3)/median(T1synthetic(grayMatter_1))
   wmNormRaw  = T1wRawMed(4)/median(T1synthetic(whiteMatter_1))
   
   std([CSFnormRaw gmNormRaw wmNormRaw])
   
   CSFnormBf = T1wBfMed(2)/median(T1synthetic(CSF_1))
   GMnormBf  = T1wBfMed(3)/median(T1synthetic(grayMatter_1))
   WMnormBf  = T1wBfMed(4)/median(T1synthetic(whiteMatter_1))
   
   std([CSFnormBf GMnormBf WMnormBf])
   
   % iterative maps
   PD_mapi = PD_map0;   
   T1_mapi = T1_map0;
   T2star_mapi = T2star_map0;   
   PD_tissueAvei = PD_tissueAve0;
   T1_tissueAvei = T1_tissueAve0;
   T2star_tissueAvei = T2star_tissueAve0;   
   
   % write to file
   if write2file
       PD_map = PD_map1;
       T1_map = T1_map1;
       T2_map = T2_map0;
       T2star_map = T2star_map1;
       save(datafileNameMATout,'PD_map','T1_map','T2_map','T2star_map','Mz','TRge','TEge','alpha');
   else    
   
   %% update PD
   % raw data
   PD_map = T1raw./(sin(alpha)*Mz.*((1-exp(-TRge/T1_mapi))./(1-cos(alpha)*exp(-TRge/T1_mapi))).*exp(-TEge/T2star_mapi));
   PD_map(PD_map<0) = 0;    
   PD_map(isnan(PD_map)) = 0;
   PD_map(isinf(PD_map)) = 0;  
   % impose max intensity in brain to be 1.0
   PD_map(PD_map>1) = 1.0; 
   %PDbrainMax = max(max(max(PD_map(CSF_1))))
   %PD_map = PD_map/PDbrainMax;
   
   PDrawMed = [median(PD_map(background_1)),median(PD_map(CSF_1)),median(PD_map(grayMatter_1)),median(PD_map(whiteMatter_1)),median(PD_map(fat_1)),median(PD_map(muscle_1)),median(PD_map(muskin_1)),median(PD_map(skull_1)),median(PD_map(vessels_1)),median(PD_map(fat2_1)),median(PD_map(dura_1)),median(PD_map(dura_1))] 

   % bilateral filtered
   PD_map = T1bf./(sin(alpha)*Mz.*((1-exp(-TRge/T1_mapi))./(1-cos(alpha)*exp(-TRge/T1_mapi))).*exp(-TEge/T2star_mapi));
   PD_map(PD_map<0) = 0;    
   PD_map(isnan(PD_map)) = 0;
   PD_map(isinf(PD_map)) = 0;  
   % impose max intensity in brain to be 1.0
   PD_map(PD_map>1) = 1.0; 
   %PDbrainMax = max(max(max(PD_map(CSF_1))))
   %PD_map = PD_map/PDbrainMax;
   
   PDbfMed = [median(PD_map(background_1)),median(PD_map(CSF_1)),median(PD_map(grayMatter_1)),median(PD_map(whiteMatter_1)),median(PD_map(fat_1)),median(PD_map(muscle_1)),median(PD_map(muskin_1)),median(PD_map(skull_1)),median(PD_map(vessels_1)),median(PD_map(fat2_1)),median(PD_map(dura_1)),median(PD_map(dura_1))] 
   
   PD_update(:,s) = (PDrawMed + PDbfMed + 3.0*PD_tissueAvei)/5.0;
   
   %% update T1
   % raw data   
   T1_map = -TRge/log((T1raw - sin(alpha)*Mz*PD_mapi.*exp(-TEge/T2star_mapi))./...
          (cos(alpha)*T1raw - sin(alpha)*Mz*PD_mapi.*exp(-TEge/T2star_mapi)));
   T1_map = real(T1_map);
   T1_map(T1_map<0) = 0;
   T1_map(T1_map>3000) = 3000.0;     
   T1_map(isnan(T1_map)) = 0;
   T1_map(isinf(T1_map)) = 0;
   
   % bilateral filtered
   T1rawMed = [median(T1_map(background_1)),median(T1_map(CSF_1)),median(T1_map(grayMatter_1)),median(T1_map(whiteMatter_1)),median(T1_map(fat_1)),median(T1_map(muscle_1)),median(T1_map(muskin_1)),median(T1_map(skull_1)),median(T1_map(vessels_1)),median(T1_map(fat2_1)),median(T1_map(dura_1)),median(T1_map(marrow_1))] 
   
   T1_map = -TRge/log((T1bf - sin(alpha)*Mz*PD_mapi.*exp(-TEge/T2star_mapi))./...
          (cos(alpha)*T1bf - sin(alpha)*Mz*PD_mapi.*exp(-TEge/T2star_mapi)));
   T1_map = real(T1_map);
   T1_map(T1_map<0) = 0;
   T1_map(T1_map>3000) = 3000.0;     
   T1_map(isnan(T1_map)) = 0;
   T1_map(isinf(T1_map)) = 0;
   
   T1bfMed = [median(T1_map(background_1)),median(T1_map(CSF_1)),median(T1_map(grayMatter_1)),median(T1_map(whiteMatter_1)),median(T1_map(fat_1)),median(T1_map(muscle_1)),median(T1_map(muskin_1)),median(T1_map(skull_1)),median(T1_map(vessels_1)),median(T1_map(fat2_1)),median(T1_map(dura_1)),median(T1_map(dura_1))] 
   
   T1_update(:,s) = (T1rawMed + T1bfMed + 3.0*T1_tissueAvei)/5.0;
   
   %% update T2*
   % raw data   
   T2star_map = real(-TEge/log((T1raw.*(1 - cos(alpha).*exp(-TRge./T1_mapi)))./...
                       (sin(alpha)*Mz*PD_mapi.*(1- exp(-TRge./T1_mapi)))));
   T2star_map(T2star_map < 0) = 0;
   T2star_map(T2star_map > T2_map0) = T2_map0(T2star_map > T2_map0);    
   T2star_map(isnan(T2star_map)) = 0;
   T2star_map(isinf(T2star_map)) = 0;
   
   T2starRawMed = [mean(T2star_map(background_1)),mean(T2star_map(CSF_1)),mean(T2star_map(grayMatter_1)),median(T2star_map(whiteMatter_1)),median(T2star_map(fat_1)),mean(T2star_map(muscle_1)),mean(T2star_map(muskin_1)),mean(T2star_map(skull_1)),mean(T2star_map(vessels_1)),mean(T2star_map(fat2_1)),mean(T2star_map(dura_1)),mean(T2star_map(marrow_1))] 
   
   % bilateral filtered
   T2star_map = real(-TEge/log((T1bf.*(1 - cos(alpha).*exp(-TRge./T1_mapi)))./...
                       (sin(alpha)*Mz*PD_mapi.*(1- exp(-TRge./T1_mapi)))));
   T2star_map(T2star_map < 0) = 0;
   T2star_map(T2star_map > T2_map0) = T2_map0(T2star_map > T2_map0);    
   T2star_map(isnan(T2star_map)) = 0;
   T2star_map(isinf(T2star_map)) = 0;
   
   T2starBfMed = [mean(T2star_map(background_1)),mean(T2star_map(CSF_1)),mean(T2star_map(grayMatter_1)),median(T2star_map(whiteMatter_1)),mean(T2star_map(fat_1)),mean(T2star_map(muscle_1)),mean(T2star_map(muskin_1)),mean(T2star_map(skull_1)),mean(T2star_map(vessels_1)),mean(T2star_map(fat2_1)),mean(T2star_map(dura_1)),mean(T2star_map(dura_1))] 
      
   T2star_update(:,s) = (T2starRawMed + T2starBfMed + 3.0*T2star_tissueAvei)/5.0;
   end
end

if ~write2file
PD_optimized = sum(PD_update,2)/20.0;
% first run
% [0.6668    1.0000    0.8390    0.7216    0.8706    0.9996    1.0000
% 0.1919    0.7098    0.7687    1.0000    1.0000]
T1_optimized = sum(T1_update,2)/20.0;
% first run
%1.0e+03 * [0.0000    2.2732    0.8667    0.5496    0.4329    0.8439    1.9508
% 2.8600    1.3113    0.5057    2.3032    1.7669]
T2star_optimized = sum(T2star_update,2)/20.0; 
% first run
%[0.1346   65.6918   60.0761   44.9873   35.5819   37.7074   20.9447
%47.8776   27.7594   57.1185   81.3867   58.2742]
end

mean(MSE0)
mean(MSE1)
mean(MSE20)
mean(MSE21)
