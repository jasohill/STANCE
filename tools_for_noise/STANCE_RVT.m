function [RVT,FRC,TIV] = STANCE_RVT(dt,Nt,RTI,RTI_sigma,weight,A_sigma,SighTI,SighTI_sigma,impulse_seed,dz_seed)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%
% Generates a respiratory volume time-series (RVT) characterized by 
% average respiratory time RTI with standard deviation RTI_sigma for a person of 
% a specific weight (optional)
% Outputs are in units of mL. 
% FRC - Functional Residual Capacity, the lung volume left after a normal exhalation
% TIV - Typical Inspiratory Volume, the lung volume after normal inhalation.
%
%   author:         Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
%   updated:        21 NOV 2016
if nargin < 10
    dz_seed = [];
end
if nargin < 9
    impulse_seed = [];
end
if nargin < 8
    SighTI_sigma = 15;
end
if isempty(SighTI_sigma)
    SighTI_sigma = 15;
end    

if nargin < 7
    SighTI = 400;
end
if isempty(SighTI)
    SighTI = 400;
end    

if nargin < 6
    A_sigma = 0.005;
end
if isempty(A_sigma)
    A_sigma = 0.005;
end

if nargin < 5
    weight = 75.0;
end
if isempty(weight)
    weight = 75.0;
end

if nargin < 4
    RTI_sigma = 0.25;
end
if isempty(RTI_sigma)
    RTI_sigma = 0.25;
end

if nargin < 3
    RTI = 4.0;
end
if isempty(RTI)
    RTI = 4.0;
end

% In regular respiration, (A. E. Lujan, J. M. Balter and R. K. Ten Haken 
% "A method for incorporating organ motion due to breathing into 3D
% dose calculations in the liver: sensitivity to variations in motion" 
% published in Med. Phys. 30 2643–9 [2003]) 
% a single respiratory cycle is represented with a functional form:
% z(phi) = A*(cos(pi*phi)^4)

% generate the respiratory motion time-series
z_ts = make_pulse_timeseries(dt,Nt,RTI,RTI_sigma,4.0,A_sigma,SighTI,SighTI_sigma,impulse_seed,dz_seed);

% Average lung volume parameters as detailed in the thesis 
% "Simulation of an Artificial Respiratory System"
% by A. Delawari & R. Doelman (2010). [Weight = 100kg]
% Weight proportions from 
% https://en.wikipedia.org/wiki/Lung_volumes
% https://en.wikipedia.org/wiki/Lung_volumes#/media/File:Lungvolumes_Updated.png
% -----------------------------------------------------------------------
% Model lungs as an ellipsoid with dimensions a = 11, b = 11, c = 6 for
% after regular exhale V = ERV + RV ~ 3.0 L.
% NAME                                    VOLUME:[L]@100kg  [mL/kg] DeltaR
% Residual Volume (RV)                          1.5         15     -1.75
% Expirational Reserve Volume (ERV)             1.5         15
% Functional Residual Capacity (FRC) = ERV + RV 3.0         30      0.00      
% Tidal volume (typical breath) (TV or V_T)     0.7         7 
% Typical inhalation volume = V_T + ERV + RV    3.7         37     +0.58
% Inspirational Reserve Volume (IRV)            4.3         43
% Vital Capacity (VC)  = IRV + TV + ERV         6.5         65
% Total Lung Capacity (TLC) = IRV+TV+ ERV+ RV = 8.0         80     +3.3
% NOTE: average weight is 75 kg -> TLC = 6.0 L

a = 11;
b = 11;
c = 6;
deltaR = 0.58*z_ts.Data;

FRC = (4*pi*weight/300.0)*a*b*c;
TIV = (4*pi*weight/300.0)*(a+0.58)*(b+0.58)*(c+0.58);
RVT_Data = (4*pi*weight/300.0).*(a+deltaR).*(b+deltaR).*(c+deltaR);

ts_time = squeeze(dt:dt:(Nt*dt));
RVT = timeseries(squeeze(RVT_Data),ts_time);
RVT.TimeInfo.Units = 'seconds';
RVT = setuniformtime(RVT,'StartTime',dt);
RVT = setuniformtime(RVT,'Interval',dt);

end

