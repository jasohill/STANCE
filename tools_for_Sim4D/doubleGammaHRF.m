function out = doubleGammaHRF(t, varargin)
% Specifies a double-gamma variate hemodynamic response function for the 
% given time vector and parameters.

%   translators/authors: Xiangyu Liu & Dr. Jason E. Hill (post-doc fellow with CNT at TTNI)
%   date:                21 NOV 2015
%   updated:             24 MAR 2017
%   inspiration:         somewhat by neuRosim's canonicalHRF() R function [Marijke Welvaert <marijke.welvaert@gmail.com>]
%   defaults:            based on the cannocical HRF found in SPM
% 
% Arguments
% t: Time vector in seconds.
% param: List of parameters of the hemodynamic response function. The list should
%        contain the following:
%        alpha1 Delay of response relative to onset   (default: 6)
%        alpha2 Delay of undershoot relative to onset (default:16)
%        beta1 Dispersion of response                 (default:1)
%        beta2 Dispersion of undershoot               (default:1)
%        c Scale of undershoot                        (default:1/6)
% verbose If true, warnings are displayed.

%% handle variable argument options
nVarargs = length(varargin);

% set initial values
param = [];
verbose = true;

if nVarargs > 0
    if mod(nVarargs,2) >0
        error('optional arguments must come in pairs');
    end
    for i = 1:2:nVarargs
        switch varargin{i}
            case 'param'
                param = varargin{i+1};
            case 'verbose'
                verbose = varargin{i+1};
            otherwise
                error('wrong input option');
        end
    end
end

% error checking
if isempty(param) 
    if verbose == true
        warning('Default parameters for HRF are used');
    end
    alpha1 = 6;
    alpha2 = 16;
    beta1 = 1.0;
    beta2 = 1.0;
    c = 1/6.0;
else
    if ~isstruct(param)
        alpha1 = param(1);
        alpha2 = param(2);
        beta1  = param(3);
        beta2  = param(4);
        c      = param(5);    
    else
        alpha1 = param.alpha1;
        alpha2 = param.alpha2;
        beta1  = param.beta1;
        beta2  = param.beta2;
        c      = param.c;          
    end
end

alpha1_map_flag = false;
alpha2_map_flag = false;
beta1_map_flag  = false;
beta2_map_flag  = false;
c_map_flag      = false;
if isfield(alpha1,'map')
    [~,Y_alpha1] = STANCE_load_volume(alpha1.map);    
    alpha1 = alpha1.max*double(Y_alpha1)/255.0;
    alpha1_map_flag = true;
end
if isfield(alpha2,'map')
    [~,Y_alpha2] = STANCE_load_volume(alpha2.map);    
    alpha2 = alpha2.max*double(Y_alpha2)/255.0;
    alpha2_map_flag = true;
end
if isfield(beta1,'map')
    [~,Y_beta1] = STANCE_load_volume(beta1.map);
    beta1 = beta1.max*double(Y_beta1)/255.0;
    beta1_map_flag = true;
end
if isfield(beta2,'map')
    [~,Y_beta2] = STANCE_load_volume(beta2.map);
    beta2 = beta2.max*double(Y_beta2)/255.0;
    beta2_map_flag = true;
end
if isfield(c,'map')
    [~,Y_c] = STANCE_load_volume(c.map);    
    c = c.max*double(Y_c)/255.0;
    c_map_flag = true;
end

%% return the double Gamma function (canonical) HRF
% d1 = alpha1*beta1;     % time to peak of response 
% d2 = alpha2*beta2;     % time to peak of undershoot
A = 2/(1-c);    % normalize integral to 2.0 (SPM default) 

if (alpha1_map_flag) || (alpha2_map_flag) || (beta1_map_flag) || (beta2_map_flag) || (c_map_flag)
   % assume 3D parameter maps
   for i = 1:length(t)
      out(:,:,:,i) = A.*(((beta1.^alpha1).*(t(i).^(alpha1-1)).*exp(-beta1*t(i))./gamma(alpha1)) ...
                   - c.*((beta2.^alpha2).*(t(i).^(alpha2-1)).*exp(-beta2*t(i))./gamma(alpha2))); %#ok<AGROW>
   end
else
     out = A.*(gampdf(t,alpha1,(1./beta1)) ...
         - c.*gampdf(t,alpha2,(1./beta2)));
% which is the same as:
%    out = A*(((t.^(alpha1-1)).*(beta1.^alpha1).*exp(-beta1*t)./gamma(alpha1)) ...
%        - c.*((t.^(alpha2-1)).*(beta2.^alpha2).*exp(-beta2*t)./gamma(alpha2)));    
end

end