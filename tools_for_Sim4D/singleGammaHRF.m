function out = singleGammaHRF(t, varargin)
% Specifies a Gamma variate hemodynamic response function for the given 
% time vector and FWHM or parameters.
% 
%   translators/authors: Xiangyu Liu & Dr. Jason E. Hill (post-doc fellow with CNT at TTNI)
%   date:                21 NOV 2015
%   updated:             24 MAR 2017
%   inspiration:         somewhat by neuRosim's gammaHRF() R function [Marijke Welvaert <marijke.welvaert@gmail.com>]
% 
% Arguments:
% t: Time vector in seconds.
% FWHM: Full Width Half Maximum of the Gamma variate function.
% (optional) k: the scaling parameter (phase-delay)
% (optional) onset: [s] the onset delay
% (optional) theta: [s] the rate parameter that determines the time-scaling
% verbose: If true, warnings are displayed. ?
% 
% NOTE: the peak delay is given by k*theta;
%       the dispersion by pd2

%% handle variable argument options
nVarargs = length(varargin);

% set initial values
FWHM  = 0;
k     = 4;
onset = 0;
theta = [];
param = [];
verbose = true;

if nVarargs > 0
    if mod(nVarargs,2) >0
        error('optional arguments must come in pairs');
    end
    for i = 1:2:nVarargs
        switch varargin{i}
            case 'FWHM'
                FWHM = double(varargin{i+1});
            case 'k'
                k = double(varargin{i+1}); 
            case 'onset'
                onset = double(varargin{i+1});                 
            case 'theta'
                theta = double(varargin{i+1});                  
            case 'param'
                param = varargin{i+1};                
            case 'verbose'
                verbose = varargin{i+1};
            otherwise
                error('unsupported or unrecognized option');
        end
    end
end

if isempty(param) 
    if verbose == true
        warning('Some default parameters for HRF may be used');
    end
else
    if ~isstruct(param)
        FWHM = param(1);
        if length(param) > 1
        k    = param(2);
        end
        if length(param) > 2
        onset  = param(3);
        end         
        if length(param) > 3
        theta  = param(4);
        end    
    else
        if isfield(param,'FWHM')
        FWHM = param.FWHM;
        end
        if isfield(param,'k')
        k = param.k;
        end
        if isfield(param,'onset')
        onset = param.onset;
        end        
        if isfield(param,'theta')
        theta = param.theta;
        end        
        if isfield(param,'alpha')
        k = param.alpha;
        end   
        if isfield(param,'beta')
        theta = 1/param.beta;
        end           
    end
end


% error checking
if FWHM == 0
    if verbose == true
        warning('Default parameters for HRF are used');
    end
    FWHM = 4;
end

%% return the Gamma function HRF
if isempty(theta)
	theta = FWHM/gamFWHM(k);
else
	% theta predefined
end

theta_map_flag  = false;
k_map_flag  = false;
onset_map_flag      = false;
if isfield(theta,'map')
    [~,Y_theta] = STANCE_load_volume(theta.map);    
    theta = theta.max*double(Y_theta)/255.0;
    theta_map_flag = true;
end
if isfield(k,'map')
    [~,Y_k] = STANCE_load_volume(k.map);    
    k = k.max*double(Y_k)/255.0;
    k_map_flag = true;
end
if isfield(onset,'map')
    [~,Y_onset] = STANCE_load_volume(onset.map);
    onset = onset.max*double(Y_onset)/255.0;
    onset_map_flag = true;
end

if (theta_map_flag) || (k_map_flag) || (onset_map_flag)
    % assume 3D parameter maps
    for i = 1:length(t)
       out(:,:,:,i) = (1./((theta.^k).*gamma(k))).*((t(i)-onset).^(k-1)).*exp(-(t(i)-onset)./theta); %#ok<AGROW>
    end

else  
    out = gampdf(t-onset,k,theta);
end
    
end