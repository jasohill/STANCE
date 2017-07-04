function [physio_4D,RP_ts,CP_ts] = STANCE_physio_4D(fn_tissues,Nt,scan,physio,veto_factor)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%
% A signal equation solver for EPI sequences.
%
%    fn_tissue:     filename of tissue priors for target space
%    scan:          MR scan protocol structure
%    physio:        physiological noise info structure
%
% author: Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
% STANCE_physio_4D.m      updated     25 MAR 2017
%
% Currently simulated main sources of respiratory and cardiac noise.
% NOTE: Sighing. deep-breathing and breath-hold will be added in future development.
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

% handle variable inputs and default values
if nargin < 5
    veto_factor = 3; % number of sigma out to despike data
end
if nargin < 4
    physio = [];
end
if nargin < 4 || isempty(physio)
    physio.weight                  = 75.0;  % [kg] weight
    physio.lambdas                 = [0.009,  0.006, 0.02,  0.05]; %  [GM WM CSF BV] lambda values for major tissue types according to the tSNR model
    physio.Rsquared                = [[0.046, 0.022, 0.042, 0.013, 0.00325, 0.022];... % CSF Row: [RP RR CP CR BP InterCRP]
                                     [0.04,   0.021, 0.03,  0.012, 0.0025,  0.017];... % GM Row
                                     [0.042,  0.020, 0.032, 0.011, 0.00275, 0.020];... % WM Row 
                                     [0.037,  0.017, 0.079, 0.009, 0.00225, 0.040];... % BS Row
                                     [0.033,  0.02,  0.25,  0.01,  0.0025,  0.02]];    % BV Row                                     
    physio.respiratory.TI          = 4.0;   % [s]  average respiratory time interval
    physio.respiratory.sigma       = 0.25;  % [s]  the standard deviation of " 
    physio.respiratory.A_z         = 1.0;   % [cm] average chest motion height  
    physio.respiratory.A_z_sigma   = 0.005; % [cm] the standard deviation of chest motion height 
    physio.respiratory.sigh.flag   = false; % do not simulate sighing
    physio.respiratory.sigh.TI     = 400.0; % [s]  average sigh time interval 
    physio.respiratory.sigh.sigma  = 15.0;  % [s]  the standard deviation of "
    physio.respiratory.seed.impulse= [];    % the random number generator seed for respiratory impulse
    physio.respiratory.seed.dz     = [];    % the random number generator seed for chest motion   
    physio.respiratory.RRF.param   = [];    % struct/list of subject's respiratory response parameters
    physio.respiratory.lags        = [];    % spatial map of time-lags 
    physio.cardiac.TI              = 1.05;  % [s]  heart beat time interval
    physio.cardiac.IPFM.freqs      = [0.02,0.1,NaN];    % [Hz] frequencies for Integral Pulse Frequency Modulation Model (IPFM)
    physio.cardiac.IPFM.sigmas     = [0.2,0.2,NaN];     % the standard deviations of rates in terms of fractional value of (1/f) for IPFM
    physio.cardiac.IPFM.amplitudes = [1.0, 1.0, (2/3)]; % the sinusoid amplitudes for IPFM
    physio.cardiac.IPFM.seeds      = [[[],[]], [[],[]], NaN]; % the random number generator seeds for IPFM
    physio.cardiac.CRF.param       = [];    % struct/list of subject's cardiac response parameters
    physio.cardiac.PWV.width       = 0.5;   % [fraction of TI] the width of the pulse wave
    physio.cardiac.PWV.width_sigma = 0.02;  % [fraction of TI] the standard deviation of the width of the pulse wave
    physio.cardiac.PWV.v0          = 3.0;   % [mm/2] the average velocity of the pulse wave in the capillaries
    physio.cardiac.PWV.v_constant_fraction = 0.1; % the fraction of the total velocity that is always flowing (minimum blood flow)
    physio.cardiac.PWV.seed        = [];  % the random number generator seed for the PWV simulator
    physio.cardiac.lags            = [];    % spatial map of time-lags
    physio.max_lag_time            = 5;     % [s] maximum lag time expected 
end
if ~isfield(physio,'weight')
    physio.weight                  = 75.0;  % [kg] weight
end
if ~isfield(physio,'lambdas')
    physio.lambdas                 = [0.009,0.006,0.02,0.05]; % lambda values for major tissue types according to the tSNR model
end
if ~isfield(physio,'Rsquared')
    physio.Rsquared                = [[0.046, 0.022, 0.042, 0.013, 0.00325, 0.022];... % CSF Row: [RP RR CP CR BP InterCRP]
                                     [0.04,   0.021, 0.03,  0.012, 0.0025,  0.017];... % GM Row
                                     [0.042,  0.020, 0.032, 0.011, 0.00275, 0.020];... % WM Row 
                                     [0.037,  0.017, 0.079, 0.009, 0.00225, 0.040];... % BS Row
                                     [0.033,  0.02,  0.25,  0.01,  0.0025,  0.02]];    % BV Row                                     
end
if ~isfield(physio,'respiratory')
    physio.respiratory.TI          = 4.0;   % [s]  average respiratory time interval
    physio.respiratory.sigma       = 0.25;  % [s]  the standard deviation of " 
    physio.respiratory.A_z         = 1.0;   % [cm] average chest motion height  
    physio.respiratory.A_z_sigma   = 0.005; % [cm] the standard deviation of chest motion height 
    physio.respiratory.sigh.flag   = false; % do not simulate sighing
    physio.respiratory.sigh.TI     = 400.0; % [s]  average sigh time interval 
    physio.respiratory.sigh.sigma  = 15.0;  % [s]  the standard deviation of "
    physio.respiratory.seed.impulse= [];    % the random number generator seed for respiratory impulse
    physio.respiratory.seed.dz     = [];    % the random number generator seed for chest motion   
    physio.respiratory.RRF.param   = [];    % struct/list of subject's respiratory response parameters
    physio.respiratory.lags        = [];    % spatial map of time-lags 
end
if ~isfield(physio.respiratory,'TI')
    physio.respiratory.TI          = 4.0;   % [s]  average respiratory time interval
end
if ~isfield(physio.respiratory,'sigma')
    physio.respiratory.sigma       = 0.25;  % [s]  the standard deviation of " 
end
if ~isfield(physio.respiratory,'A_z')
    physio.respiratory.A_z         = 0.1;   % [cm] average chest motion height 
end
if ~isfield(physio.respiratory,'A_z_sigma')
    physio.respiratory.A_z_sigma   = 0.005;   % [cm] the standard deviation of chest motion height 
end
if ~isfield(physio.respiratory,'RRF')
    physio.respiratory.RRF.param   = [];    % struct/list of subject's respiratory response parameters
end
if ~isfield(physio.respiratory,'sigh')
    physio.respiratory.sigh.flag  = false; % do not simulate sighing
    physio.respiratory.sigh.TI    = 400.0; % [s]  average sigh time interval 
    physio.respiratory.sigh.sigma = 15.0;  % [s]  the standard deviation of "
end
if ~isfield(physio.respiratory.sigh,'flag')
    physio.respiratory.sigh.flag   = false; % do not simulate sighing
end
if ~isfield(physio.respiratory.sigh,'TI')
    physio.respiratory.sigh.TI     = 400.0; % [s]  average sigh time interval 
end
if ~isfield(physio.respiratory.sigh,'sigma')
    physio.respiratory.sigh.sigma = 15.0;  % [s]  the standard deviation of "
end
if ~isfield(physio.respiratory,'lags')
    physio.respiratory.lags       = [];  % [integer increments of dt] time-series lag times
end
if ~isfield(physio.respiratory,'seed')
    physio.respiratory.seed.impulse= [];    % the random number generator seed for respiratory impulse
    physio.respiratory.seed.dz     = [];    % the random number generator seed for chest motion   
end
if ~isfield(physio.respiratory.seed,'impulse')
    physio.respiratory.seed.impulse= [];    % the random number generator seed for respiratory impulse
end
if ~isfield(physio.respiratory.seed,'dz')
    physio.respiratory.seed.dz     = [];    % the random number generator seed for chest motion   
end
if ~isfield(physio,'cardiac')
    physio.cardiac.TI              = 1.05;  % [s]  heart beat time interval
    physio.cardiac.IPFM.freqs      = [0.02,0.1,NaN];    % [Hz] frequencies for Integral Pulse Frequency Modulation Model (IPFM)
    physio.cardiac.IPFM.sigmas     = [0.2,0.2,NaN];     % the standard deviation of rates in terms of fractional value of (1/f) for IPFM
    physio.cardiac.IPFM.amplitudes = [1.0, 1.0, (2/3)]; % the standard deviation of sinusoid amplitudes  for IPFM
    physio.cardiac.IPFM.seeds      = [[[],[]], [[],[]], NaN]; % the random number generator seeds for IPFM
    physio.cardiac.PWV.width       = 0.5;   % [fraction of TI] the width of the pulse wave
    physio.cardiac.PWV.width_sigma = 0.02;  % [fraction of TI] the standard deviation of the width of the pulse wave
    physio.cardiac.PWV.v0          = 3.0;   % [mm/2] the average velocity of the pulse wave in the capillaries
    physio.cardiac.PWV.v_constant_fraction = 0.1; % the fraction of the total velocity that is always flowing (minimum blood flow)
    physio.cardiac.PWV.seed        = [];  % the random number generator seed for the PWV simulator
    physio.cardiac.CRF.param       = [];    % struct/list of subject's cardiac response parameters
end
if ~isfield(physio.cardiac,'TI')
    physio.cardiac.TI              = 1.05;              % [s]  heart beat time intervalend
end
if ~isfield(physio.cardiac,'IPFM')
    physio.cardiac.IPFM.freqs      = [0.02,0.1,NaN];    % [Hz] frequencies for Integral Pulse Frequency Modulation Model (IPFM)
    physio.cardiac.IPFM.sigmas     = [0.2,0.2,NaN];     % the standard deviation of rates in terms of fractional value of (1/f) for IPFM
    physio.cardiac.IPFM.amplitudes = [1.0, 1.0, (2/3)]; % the standard deviation of sinusoid amplitudes  for IPFM
    physio.cardiac.IPFM.seeds      = [[[],[]], [[],[]], NaN]; % the random number generator seeds for IPFM
end
if ~isfield(physio.cardiac.IPFM,'freqs')
    physio.cardiac.IPFM.freqs      = [0.02,0.1,NaN];    % [Hz] frequencies for Integral Pulse Frequency Modulation Model (IPFM)
end
if ~isfield(physio.cardiac.IPFM,'sigmas')
    physio.cardiac.IPFM.sigmas     = [0.2,0.2,NaN];     % the standard deviation of rates in terms of fractional value of (1/f) for IPFM
end
if ~isfield(physio.cardiac.IPFM,'amplitudes')
    physio.cardiac.IPFM.amplitudes = [1.0, 1.0, (2/3)]; % the standard deviation of sinusoid amplitudes  for IPFM
end
if ~isfield(physio.cardiac.IPFM,'seeds')
    physio.cardiac.IPFM.seeds      = [[[],[]], [[],[]], NaN]; % the random number generator seeds for IPFM
end
if ~isfield(physio.cardiac,'PWV')
    physio.cardiac.PWV.width       = 0.5;   % [fraction of TI] the width of the pulse wave
    physio.cardiac.PWV.width_sigma = 0.02;  % [fraction of TI] the standard deviation of the width of the pulse wave
    physio.cardiac.PWV.v0          = 3.0;   % [mm/2] the average velocity of the pulse wave in the capillaries
    physio.cardiac.PWV.v_constant_fraction = 0.1; % the fraction of the total velocity that is always flowing (minimum blood flow)
    physio.cardiac.PWV.seed        = [];  % the random number generator seed for the PWV simulator
end
if ~isfield(physio.cardiac.PWV,'width')
    physio.cardiac.PWV.width       = 0.5;   % [fraction of TI] the width of the pulse wave
end
if ~isfield(physio.cardiac.PWV,'width_sigma')
    physio.cardiac.PWV.width_sigma = 0.02;  % [fraction of TI] the standard deviation of the width of the pulse wave
end
if ~isfield(physio.cardiac.PWV,'v0')
    physio.cardiac.PWV.v0          = 3.0;   % [mm/2] the average velocity of the pulse wave in the capillaries
end
if ~isfield(physio.cardiac.PWV,'v_constant_fraction')
    physio.cardiac.PWV.v_constant_fraction = 0.1; % the fraction of the total velocity that is always flowing (minimum blood flow)
end
if ~isfield(physio.cardiac.PWV,'seed')
    physio.cardiac.PWV.seed        = [];  % the random number generator seed for the PWV simulator
end
if ~isfield(physio.cardiac,'CRF')
    physio.cardiac.CRF.param       = [];    % struct/list of subject's cardiac response parameters
end
if ~isfield(physio.cardiac,'lags')
    physio.cardiac.lags = [];  % [integer increments of dt] time-series lag times
end
if ~isfield(physio,'max_lag_time')
    physio.max_lag_time = 5;  % [s] maximum lag time (+/-) expected
end
% For a reference on the tSNR model of fMRI noise see "Physiological noise in oxygenation-sensitive magnetic resonance imaging" 
% by G. Krüger & G. H. Glover in Magn Reson Med. 2001 Oct;46(4):631-7.
% Rsquared values from Fig. 2 in "Hemodynamic response function in resting
% brain:  disambiguating neural events and autonomic effects"
% by Guo Guo-Rong & Daniele Marinazzo available at bioR?iv (beta): the preprint server for biology, doi:
% https://doi.org/10.1101/028514 (2015).
% NOTE: estimate for blood vessel Rsquared comes from average values plus 
% the fact that the blood flow velocity in a blood vessel is ~70 times that 
% of capillaries in the GM so the CP std is expected to be ~8 times larger.
% NOTE: Mayer wave Rsquared estimated at ~25% of that due to heart rate variability (CR) 
% see Katura, T., Tanaka, N., Obata, A., Sato, H., and Makia. A., 
% "Quantitative evaluation of interrelations between spontaneous low frequency 
% oscillations in cerebral hemodynamics and systemic cardiovascular dynamics," 
% NeuroImage 31(4), 1592-1600 (2006).

if nargin < 3
    scan = [];
end
if nargin < 3 || isempty(scan)
    scan.voxel.size    = [3 3 3]; % [mm] voxel size dimensions
    scan.voxel.matrix  = [64 64 NaN]; % size of matrix (number of Z slices may be autodetected)
    scan.voxel.spacing = [0 0 0.2*scan.voxel.size(3)]; % assume 20% Z spacing
    scan.tiltAngle     = 15;   % [degrees] tilt angle from AC-PC line
    scan.TR            = 2000; % [ms] repetition time
    scan.TE            = 30;   % [ms] echo time
    scan.ES            = 0.51; % [ms] echo spacing
    scan.FA            = 90;   % [degrees] flip angle 
    scan.BW            = 2232; % [Hz/Px] bandwidth
    scan.order         = 'SD'; % SD = sequential descending order
    scan.B0             = 3.0; % [T] the main magnet field strength    
    scan.KM0           = 2225; % fit to data with max of 909 at 3T and FA = 90 degrees
    scan.noise_method  = 'sigma';
    scan.noise         = 0;   
    scan.attenuation   = 0;
    scan.acceleration  = 1;
end
if ~isfield(scan,'voxel')
    scan.voxel.size    = [3 3 3]; % [mm] voxel size dimensions
    scan.voxel.matrix  = [64 64 NaN];  % size of matrix (number of Z slices may be autodetected)
    scan.voxel.spacing = [0 0 0.2*scan.voxel.size(3)]; % assume 20% Z spacing
end
if ~isfield(scan.voxel,'size')
    scan.voxel.size    = [3 3 3]; % [mm] voxel size dimensions
end
if ~isfield(scan.voxel,'matrix')
    scan.voxel.matrix  = [64 64 NaN]; % size of matrix (number of Z slices may be autodetected)
end
if ~isfield(scan.voxel,'spacing')
    scan.voxel.spacing = [0 0 0.2*scan.voxel.size(3)]; % assume 20% Z spacing
end
if ~isfield(scan,'tiltAngle')
    scan.tiltAngle     = 15; % [degrees] 
end
if ~isfield(scan,'TR')
    scan.TR            = 2000; % [ms]
end
if ~isfield(scan,'TE')
    scan.TE            = 30;   % [ms]
end
if ~isfield(scan,'ES')
    scan.ES            = 0.51; % [ms] echo spacing
end
if ~isfield(scan,'FA')
    scan.FA            = 90;   % degrees 
end
if ~isfield(scan,'BW') 
    scan.BW            = 2232; % [Hz/Px] bandwidth
end
if ~isfield(scan,'order') 
    scan.order         = 'SD'; % SD = sequential descending order
end
if ~isfield(scan,'B0') 
    scan.B0             = 3.0; % [T] the main magnet field strength    
end
if ~isfield(scan,'KM0') 
    scan.KM0           = 2225; % fit to data with max of 909 at 3T and FA = 90 degrees
end
if ~isfield(scan,'noise_method') 
    scan.noise_method  = 'sigma';
end
if ~isfield(scan,'noise_method') 
    scan.noise_method  = 'sigma';
end
if ~isfield(scan,'noise') 
    scan.noise         = 0;   
end
if ~isfield(scan,'attenuation')   
    scan.attenuation   = 0;
end
if ~isfield(scan,'acceleration') 
    scan.acceleration  = 1;
end

% define parameters describing the subjects physiological signals and the
% MR scan protocol
weight         = physio.weight;
lambdas        = physio.lambdas;
RTI            = physio.respiratory.TI;
RTI_sigma      = physio.respiratory.sigma;
A_z_sigma      = physio.respiratory.A_z_sigma;
RRF_param      = physio.respiratory.RRF.param;
STI            = physio.respiratory.sigh.TI; 
STI_sigma      = physio.respiratory.sigh.sigma;
impulse_seed   = physio.respiratory.seed.impulse;
dz_seed        = physio.respiratory.seed.dz;  
TI             = physio.cardiac.TI;
IPFM_fs        = physio.cardiac.IPFM.freqs;
IPFM_sigmas    = physio.cardiac.IPFM.sigmas;
IPFM_As        = physio.cardiac.IPFM.amplitudes;
IPFM_seeds     = physio.cardiac.IPFM.seeds;
PWV_width      = physio.cardiac.PWV.width;
PWV_width_sigma= physio.cardiac.PWV.width_sigma;
PWV_v0         = physio.cardiac.PWV.v0; 
PWV_v_constant = physio.cardiac.PWV.v_constant_fraction;
PWV_seed       = physio.cardiac.PWV.seed;       
CRF_param      = physio.cardiac.CRF.param;
lags_R         = physio.respiratory.lags;
lags_C         = physio.cardiac.lags;
if isempty(lags_R) && isempty(lags_C)
    max_lag    = 0;
else
    max_lag    = physio.max_lag_time;
end
CSF_Rsquared = physio.Rsquared(1,:);
GM_Rsquared  = physio.Rsquared(2,:);
WM_Rsquared  = physio.Rsquared(3,:);
BS_Rsquared  = physio.Rsquared(4,:);
BV_Rsquared  = physio.Rsquared(5,:);
%ave_Rsquared = mean(physio.Rsquared);

TR         = scan.TR;
TE         = scan.TE;
FA         = scan.FA;
order      = scan.order;
KM0        = scan.KM0;
method     = scan.noise_method;
noise      = scan.noise;
attenuation = scan.attenuation;


if ~exist('STANCE.mat') %#ok<*EXIST>
    if ~exist('../STANCE.mat')
       load('../../STANCE.mat');
    else
       load('../STANCE.mat');
    end
else
    load('STANCE.mat');
end

tic

[~,Y_Tissue_Labels] = STANCE_load_volume(fn_tissues);

%% design 4d timeseries

Nslices = size(Y_Tissue_Labels,3);
TRsec = TR/1000;
dt = TRsec/Nslices;
%Nt = length(exp_design.Data);  % length of fMRI experiment
Nt_SS = ceil(60/dt);           % add 60 seconds to beginning of physiological time series to allow RR and CR to reach steady state
Nt_max_lag = ceil(max_lag/dt); % pad time series with this ammount before and after
times = [dt*((Nt_SS+Nt_max_lag):-1:0) dt*(1:Nt) dt*(Nt+(1:Nt_max_lag))]';
prepadding = length(Nt_SS:-1:0);
% NV = (Nt/Nslices);     % number of volumes to be generated
Nphys = length(times); % length of physiological noise simulation


%% simulate breathing-related physiological time-series

% respiratory volume time-series
display('o Generating the respiratory volume time-series (RVT).')
[RVT,FRC,TIV] = STANCE_RVT(dt,Nphys,RTI,RTI_sigma,weight,A_z_sigma,STI,STI_sigma,impulse_seed,dz_seed);
RVT_Data = squeeze(RVT.Data);

% Define RVT baseline to use in modelling HRV
RVT_baseline_Data = (RVT_Data - (FRC+0.5*(TIV-FRC)))/(0.5*(TIV-FRC));
RVT_baseline = timeseries(squeeze(RVT_baseline_Data),times);

% Define the respiratory pulse (RP) time-series with 0 average and unit the standard deviation 
display('o Generating the respiratory pulse (RP) time-series.')
RP_ts = RVT_baseline;
RP_ts.Data = (RVT_baseline_Data - mean(RVT_baseline_Data(Nt_SS+2:end)))./std(RVT_baseline_Data(Nt_SS+2:end));
% Increased lung volume increases magnetic field inhomogeneity, which DECREASES T2*.
RP_ts_Data = -squeeze(RP_ts.Data);

% Define the respiratory response (RR) time-series with 0 average and unit the standard deviation
display('o Generating the respiratory response (RR) time-series.')
RR_ts = STANCE_apply_response_function(dt, RVT_baseline,'RRF', RRF_param);
RR_ts.Data = (RR_ts.Data - mean(RR_ts.Data(Nt_SS+2:end)))./std(RR_ts.Data(Nt_SS+2:end));
RR_ts_Data = squeeze(RR_ts.Data);


%% simulate cardiac-related physiological time-series

% This section follows the Integral Pulse Frequency Modulation (IPFM) model. 
% Model heart-rate variability with realistic quasi-periodic functions 
% from physiological sources of breathing and blood pressure variations.
% Three tone harmonic modulation is in the well-known HRV spectral regions of very low
% frequency (VLF), low frequency (LF) and high frequency (HF), e.g. see in section 4.1 of
% "Spectral analysis of heart rate variability using the integral pulse frequency modulation model"
% by I. P. Mitov, published in Med. Biol. Eng. Comput., 2001, 39, 1-7.
TI = 1.05; % beat/second
m0 = dt./TI;

% Define the very low frequency (VLF) heart-rate modulation time-series
display('o Generating the VLF heart-rate modulation time-series.')
f1 = IPFM_fs(1); %Hz
VLF = make_pulse_timeseries(dt,Nphys,1/f1,IPFM_sigmas(1)*(1/f1),2);
VLF_Data = squeeze(VLF.Data);
VLF_Data = 2.0*(VLF_Data - 0.5);

% Define the Mayer wave blood pressure variation time-series
display('o Generating Mayer wave blood pressure (BP) variation time-series.')
f2 = IPFM_fs(2); %Hz
Mayer_wave = make_pulse_timeseries(dt,Nphys,1/f2,IPFM_sigmas(2)*(1/f2),2);
Mayer_wave_Data = squeeze(Mayer_wave.Data);
Mayer_wave_Data = 2.0*(Mayer_wave_Data - 0.5);

m_phys = m0*(IPFM_As(1)*VLF_Data+IPFM_As(2)*Mayer_wave_Data+IPFM_As(3)*RVT_baseline_Data);

ms_data = m0 + m_phys;
ms = timeseries(ms_data,times);
ms.TimeInfo.Units = 'seconds';
ms = setuniformtime(ms,'StartTime',dt);
ms = setuniformtime(ms,'Interval',dt);

% Instantaneous heart rate (HR) time-series
display('o Generating the instantaneous heart rate (HR) time-series from a generalized IPFM model.')
HR_data = ms_data*TI*60./dt;
HR_ts = timeseries(squeeze(HR_data),times);

% Cardiac response (CR) time-series
display('o Generating the cardiac response (CR) time-series from HR.')
CR_ts = STANCE_apply_response_function(dt, HR_ts, 'CRF', CRF_param);
CR_ts.Data = (CR_ts.Data - mean(CR_ts.Data(Nt_SS+2:end)))./std(CR_ts.Data(Nt_SS+2:end));
CR_ts_Data = squeeze(CR_ts.Data);

% Cardiac event impulse (CEI) time-series
display('o Generating the cardiac event impulse (CEI) time-series.')
[CEI_ts,~] = make_event_timeseries(ms,TI);
CEI_ts_Data = squeeze(CEI_ts.Data);

% Pulse wave velocity (PWV) time-series
display('o Generating the pulse wave velocity (PWV) time-series.')
v_ts = STANCE_PWV_timeseries(CEI_ts,PWV_width,PWV_width_sigma,PWV_v0,PWV_v_constant,PWV_seed);
v_ts_Data = squeeze(v_ts.Data);

% Define the cardiac pulse (CP) time-series
display('o Generating the cardiac pulse (CP) time-series.')
CP_ts = v_ts;
CP_ts.Data = (CP_ts.Data - mean(CP_ts.Data(Nt_SS+2:end)))./std(CP_ts.Data(Nt_SS+2:end));
% Cardiac pulse blood flow causes the signal S_0 to DECREASE, but T2* to INCREASE via the CRF (HgO_2 pulsations)
CP_ts_Data = -squeeze(CP_ts.Data);

% Define the interaction of the cardiac and respiratory pulses (InterCRP) time-series
display('o Generating the interaction of the cardiac and respiratory pulses (InterCRP) time-series.')
% model as a product of a weighted average of respiratory sources with cardiac sources
% ... in gray matter
InterCRP_GM_ts = RVT;
InterCRP_GM_ts.Data = (GM_Rsquared(1)*RP_ts_Data+GM_Rsquared(2)*RR_ts_Data).*(GM_Rsquared(3)*CP_ts_Data + GM_Rsquared(4).*CR_ts_Data + GM_Rsquared(5)*Mayer_wave_Data);
InterCRP_GM_mean = mean(InterCRP_GM_ts.Data(Nt_SS+2:end));
InterCRP_GM_std  = std(InterCRP_GM_ts.Data(Nt_SS+2:end));
InterCRP_GM_ts.Data = (InterCRP_GM_ts.Data - InterCRP_GM_mean)./InterCRP_GM_std;
InterCRP_GM_ts_Data = squeeze(InterCRP_GM_ts.Data);
% ... in white matter
InterCRP_WM_ts = RVT;
InterCRP_WM_ts.Data = (WM_Rsquared(1)*RP_ts_Data+WM_Rsquared(2)*RR_ts_Data).*(WM_Rsquared(3)*CP_ts_Data + WM_Rsquared(4).*CR_ts_Data + WM_Rsquared(5)*Mayer_wave_Data);
InterCRP_WM_mean    = mean(InterCRP_WM_ts.Data(Nt_SS+2:end));
InterCRP_WM_std      = std(InterCRP_WM_ts.Data(Nt_SS+2:end));
InterCRP_WM_ts.Data = (InterCRP_WM_ts.Data - InterCRP_WM_mean)./InterCRP_WM_std;
InterCRP_WM_ts_Data = squeeze(InterCRP_WM_ts.Data);
% ... in cerebrospinal fluid
InterCRP_CSF_ts = RVT;
InterCRP_CSF_ts.Data = (CSF_Rsquared(1)*RP_ts_Data+CSF_Rsquared(2)*RR_ts_Data).*(CSF_Rsquared(3)*CP_ts_Data + CSF_Rsquared(4).*CR_ts_Data + CSF_Rsquared(5)*Mayer_wave_Data);
InterCRP_CSF_mean    = mean(InterCRP_CSF_ts.Data(Nt_SS+2:end));
InterCRP_CSF_std     = std(InterCRP_CSF_ts.Data(Nt_SS+2:end));
InterCRP_CSF_ts.Data = (InterCRP_CSF_ts.Data - InterCRP_CSF_mean)./InterCRP_CSF_std;
InterCRP_CSF_ts_Data = squeeze(InterCRP_CSF_ts.Data);
% ... in blood vessels
InterCRP_BV_ts = RVT;
InterCRP_BV_ts.Data = (BV_Rsquared(1)*RP_ts_Data+BV_Rsquared(2)*RR_ts_Data).*(BV_Rsquared(3)*CP_ts_Data + BV_Rsquared(4).*CR_ts_Data + BV_Rsquared(5)*Mayer_wave_Data);
InterCRP_BV_mean    = mean(InterCRP_BV_ts.Data(Nt_SS+2:end));
InterCRP_BV_std     = std(InterCRP_BV_ts.Data(Nt_SS+2:end));
InterCRP_BV_ts.Data = (InterCRP_BV_ts.Data - InterCRP_BV_mean)./InterCRP_BV_std;
InterCRP_BV_ts_Data = squeeze(InterCRP_BV_ts.Data);
% ... in brain stem
InterCRP_BS_ts = RVT;
InterCRP_BS_ts.Data = (BS_Rsquared(1)*RP_ts_Data+BS_Rsquared(2)*RR_ts_Data).*(BS_Rsquared(3)*CP_ts_Data + BS_Rsquared(4).*CR_ts_Data + BS_Rsquared(5)*Mayer_wave_Data);
InterCRP_BS_mean    = mean(InterCRP_BS_ts.Data(Nt_SS+2:end));
InterCRP_BS_std     = std(InterCRP_BS_ts.Data(Nt_SS+2:end));
InterCRP_BS_ts.Data = (InterCRP_BS_ts.Data - InterCRP_BS_mean)./InterCRP_BS_std;
InterCRP_BS_ts_Data = squeeze(InterCRP_BS_ts.Data);


%% Physiological noise per tissue type time-series

display('o Generating the physiological noise per tissue type time-series.')

% ------------------------------
% all components (no lag times)
% ------------------------------

% Gray matter
GM_phys = RVT;
GM_phys.Data = GM_Rsquared(1)*RP_ts_Data + GM_Rsquared(2)*RR_ts_Data ...
               + GM_Rsquared(3)*CP_ts_Data + GM_Rsquared(4)*CR_ts_Data ...
               + GM_Rsquared(5)*Mayer_wave_Data + GM_Rsquared(6)*InterCRP_GM_ts_Data;
GM_phys.Data = (GM_phys.Data - mean(GM_phys.Data(Nt_SS+2:end)))./std(GM_phys.Data(Nt_SS+2:end));
GM_phys_noise = squeeze(GM_phys.Data(Nt_SS+2:end));

% White matter
WM_phys = RVT;
WM_phys.Data = WM_Rsquared(1)*RP_ts_Data + WM_Rsquared(2)*RR_ts_Data ...
               + WM_Rsquared(3)*CP_ts_Data + WM_Rsquared(4)*CR_ts_Data ...
               + WM_Rsquared(5)*Mayer_wave_Data + WM_Rsquared(6)*InterCRP_WM_ts_Data;
WM_phys.Data = (WM_phys.Data - mean(WM_phys.Data(Nt_SS+2:end)))./std(WM_phys.Data(Nt_SS+2:end));
WM_phys_noise = squeeze(WM_phys.Data(Nt_SS+2:end));

% Cerebrospinal fluid
CSF_phys = RVT;
CSF_phys.Data = CSF_Rsquared(1)*RP_ts_Data + CSF_Rsquared(2)*RR_ts_Data ...
                + CSF_Rsquared(3)*CP_ts_Data + CSF_Rsquared(4)*CR_ts_Data ...
                + CSF_Rsquared(5)*Mayer_wave_Data + CSF_Rsquared(6)*InterCRP_CSF_ts_Data;
CSF_phys.Data = (CSF_phys.Data - mean(CSF_phys.Data(Nt_SS+2:end)))./std(CSF_phys.Data(Nt_SS+2:end));
CSF_phys_noise = squeeze(CSF_phys.Data(Nt_SS+2:end));

% brain stem
BS_phys = RVT;
BS_phys.Data = BS_Rsquared(1)*RP_ts_Data + BS_Rsquared(2)*RR_ts_Data ...
               + BS_Rsquared(3)*CP_ts_Data + BS_Rsquared(4)*CR_ts_Data ...
               + BS_Rsquared(5)*Mayer_wave_Data + BS_Rsquared(6)*InterCRP_BS_ts_Data;
BS_phys.Data = (BS_phys.Data - mean(BS_phys.Data(Nt_SS+2:end)))./std(BS_phys.Data(Nt_SS+2:end));
BS_phys_noise = squeeze(BS_phys.Data(Nt_SS+2:end));

% blood vessels
BV_phys = RVT;
BV_phys.Data = BV_Rsquared(1)*RP_ts_Data + BV_Rsquared(2)*RR_ts_Data ...
               + BV_Rsquared(3)*CP_ts_Data + BV_Rsquared(4)*CR_ts_Data ...
               + BV_Rsquared(5)*Mayer_wave_Data + BV_Rsquared(6)*InterCRP_BV_ts_Data;
BV_phys.Data = (BV_phys.Data - mean(BV_phys.Data(Nt_SS+2:end)))./std(BV_phys.Data(Nt_SS+2:end));
BV_phys_noise = squeeze(BV_phys.Data(Nt_SS+2:end));


% ------------------------------
% respiratory component
% ------------------------------

% Gray matter
GM_physR = RVT;
GM_physR.Data = GM_Rsquared(1)*RP_ts_Data + GM_Rsquared(2)*RR_ts_Data;
GM_physR_raw = squeeze(GM_physR.Data(Nt_SS+2:end));
GM_physR.Data = (GM_physR.Data - mean(GM_physR.Data(Nt_SS+2:end)))./std(GM_physR.Data(Nt_SS+2:end));
GM_physR_noise = squeeze(GM_physR.Data(Nt_SS+2:end));

% White matter
WM_physR = RVT;
WM_physR.Data = WM_Rsquared(1)*RP_ts_Data + WM_Rsquared(2)*RR_ts_Data;
WM_physR_raw = squeeze(WM_physR.Data(Nt_SS+2:end));
WM_physR.Data = (WM_physR.Data - mean(WM_physR.Data(Nt_SS+2:end)))./std(WM_physR.Data(Nt_SS+2:end));
WM_physR_noise = squeeze(WM_physR.Data(Nt_SS+2:end));

% Cerebrospinal fluid
CSF_physR = RVT;
CSF_physR.Data = CSF_Rsquared(1)*RP_ts_Data + CSF_Rsquared(2)*RR_ts_Data;
CSF_physR_raw = squeeze(CSF_physR.Data(Nt_SS+2:end));
CSF_physR.Data = (CSF_physR.Data - mean(CSF_physR.Data(Nt_SS+2:end)))./std(CSF_physR.Data(Nt_SS+2:end));
CSF_physR_noise = squeeze(CSF_physR.Data(Nt_SS+2:end));

% brain stem
BS_physR = RVT;
BS_physR.Data = BS_Rsquared(1)*RP_ts_Data + BS_Rsquared(2)*RR_ts_Data;
BS_physR_raw = squeeze(BS_physR.Data(Nt_SS+2:end));
BS_physR.Data = (BS_physR.Data - mean(BS_physR.Data(Nt_SS+2:end)))./std(BS_physR.Data(Nt_SS+2:end));
BS_physR_noise = squeeze(BS_physR.Data(Nt_SS+2:end));

% blood vessels
BV_physR = RVT;
BV_physR.Data = BV_Rsquared(1)*RP_ts_Data + BV_Rsquared(2)*RR_ts_Data;
BV_physR_raw = squeeze(BV_physR.Data(Nt_SS+2:end));
BV_physR.Data = (BV_physR.Data - mean(BV_physR.Data(Nt_SS+2:end)))./std(BV_physR.Data(Nt_SS+2:end));
BV_physR_noise = squeeze(BV_physR.Data(Nt_SS+2:end));


% ------------------------------
% cardiac component
% ------------------------------

% Gray matter
GM_physC = RVT;
GM_physC.Data = GM_Rsquared(3)*CP_ts_Data + GM_Rsquared(4)*CR_ts_Data + GM_Rsquared(5)*Mayer_wave_Data;
GM_physC_raw = squeeze(GM_physC.Data(Nt_SS+2:end));
GM_physC.Data = (GM_physC.Data - mean(GM_physC.Data(Nt_SS+2:end)))./std(GM_physC.Data(Nt_SS+2:end));
GM_physC_noise = squeeze(GM_physC.Data(Nt_SS+2:end));

% White matter
WM_physC = RVT;
WM_physC.Data = WM_Rsquared(3)*CP_ts_Data + WM_Rsquared(4)*CR_ts_Data + WM_Rsquared(5)*Mayer_wave_Data;
WM_physC_raw = squeeze(WM_physC.Data(Nt_SS+2:end));
WM_physC.Data = (WM_physC.Data - mean(WM_physC.Data(Nt_SS+2:end)))./std(WM_physC.Data(Nt_SS+2:end));
WM_physC_noise = squeeze(WM_physC.Data(Nt_SS+2:end));

% Cerebrospinal fluid
CSF_physC = RVT;
CSF_physC.Data = CSF_Rsquared(3)*CP_ts_Data + CSF_Rsquared(4)*CR_ts_Data + CSF_Rsquared(5)*Mayer_wave_Data;
CSF_physC_raw = squeeze(CSF_physC.Data(Nt_SS+2:end));
CSF_physC.Data = (CSF_physC.Data - mean(CSF_physC.Data(Nt_SS+2:end)))./std(CSF_physC.Data(Nt_SS+2:end));
CSF_physC_noise = squeeze(CSF_physC.Data(Nt_SS+2:end));

% brain stem
BS_physC = RVT;
BS_physC.Data = BS_Rsquared(3)*CP_ts_Data + BS_Rsquared(4)*CR_ts_Data + BS_Rsquared(5)*Mayer_wave_Data;
BS_physC_raw = squeeze(BS_physC.Data(Nt_SS+2:end));
BS_physC.Data = (BS_physC.Data - mean(BS_physC.Data(Nt_SS+2:end)))./std(BS_physC.Data(Nt_SS+2:end));
BS_physC_noise = squeeze(BS_physC.Data(Nt_SS+2:end));

% blood vessels
BV_physC = RVT;
BV_physC.Data = BV_Rsquared(3)*CP_ts_Data + BV_Rsquared(4)*CR_ts_Data + BV_Rsquared(5)*Mayer_wave_Data;
BV_physC_raw = squeeze(BV_physC.Data(Nt_SS+2:end));
BV_physC.Data = (BV_physC.Data - mean(BV_physC.Data(Nt_SS+2:end)))./std(BV_physC.Data(Nt_SS+2:end));
BV_physC_noise = squeeze(BV_physC.Data(Nt_SS+2:end));


%% Define fuzzy memberships
%---------------------------

background  = Y_Tissue_Labels(:,:,:,1);
CSF         = Y_Tissue_Labels(:,:,:,2);
grayMatter  = Y_Tissue_Labels(:,:,:,3);
whiteMatter = Y_Tissue_Labels(:,:,:,4);
vessels     = Y_Tissue_Labels(:,:,:,9);

sliceOrder = scan.order;
sliceTiming = make_slice_timing(sliceOrder,Nslices);


hw = waitbar(0,'Computing physiological noise timeseries ...');
if isempty(lags_R) && isempty(lags_C)
    for t = (1:Nt)
        ti = ceil(t/Nslices);
        STi = mod(t,Nslices);
        if STi == 0
            STi = Nslices;
        end
        k = sliceTiming(STi);
        physio_4D(:,:,k,ti) = lambdas(1)*grayMatter(:,:,k)*GM_phys_noise(t) ...
            + lambdas(2)*whiteMatter(:,:,k)*WM_phys_noise(t) ...
            + lambdas(3)*CSF(:,:,k)*CSF_phys_noise(t) ...
            + lambdas(4)*vessels(:,:,k)*BV_phys_noise(t); %#ok<AGROW>   
        waitbar(t/Nt,hw)
    end
else
    display('o Adding lag times ...')
    if isempty(lags_R)
        lags_R = 0*lags_C;
    end
    if isempty(lags_C)
        lags_C = 0*lags_R;
    end    
    for t = (1:Nt) + Nt_max_lag
        ti = ceil(t/Nslices);
        STi = mod(t,Nslices);
        if STi == 0
            STi = Nslices;
        end
        k = sliceTiming(STi);
        X = size(CSF,1);
        Y = size(CSF,2);        
        for i = 1:X
            for j = 1:Y
                t_R = t - lags_R(i,j,k);            
                t_C = t - lags_C(i,j,k);             
                physio_4D(i,j,k,ti) = ...
                    lambdas(1)*grayMatter(i,j,k)*(GM_physR_noise(t_R) + GM_physC_noise(t_C) - (GM_physR_raw(t_R)*GM_physC_raw(t_C) + InterCRP_GM_mean)/InterCRP_GM_std) ...
                    + lambdas(2)*whiteMatter(i,j,k)*(WM_physR_noise(t_R) + WM_physC_noise(t_C)  - (WM_physR_raw(t_R)*WM_physC_raw(t_C) + InterCRP_WM_mean)/InterCRP_WM_std)  ...
                    + lambdas(3)*CSF(i,j,k)*(CSF_physR_noise(t_R) + CSF_physC_noise(t_C) - (CSF_physR_raw(t_R)*CSF_physC_raw(t_C) + InterCRP_CSF_mean)/InterCRP_CSF_std) ...
                    + lambdas(4)*vessels(i,j,k)*(BV_physR_noise(t_R) + BV_physC_noise(t_C) - (BV_physR_raw(t_R)*BV_physC_raw(t_C) + InterCRP_BV_mean)/InterCRP_BV_std);                      %#ok<AGROW>
            end
        end
       waitbar(t/Nt,hw)        
    end  
end
% despike
physio_sigma = max(lambdas);

physio_4D(physio_4D>veto_factor*physio_sigma) = veto_factor*physio_sigma;
physio_4D(physio_4D<-veto_factor*physio_sigma) = -veto_factor*physio_sigma;

close(hw)

cd(STANCEroot);



end

