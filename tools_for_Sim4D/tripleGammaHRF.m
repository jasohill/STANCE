function out = tripleGammaHRF(t, varargin)
% Specifies a triple-Gamma variate hemodynamic response function for the 
% given time vector and parameters.

%   author:        Dr. Jason E. Hill (post-doc fellow with CNT at TTNI)
%   updated:       24 MAR 2017
%
% Arguments
% t: Time vector in seconds.
% param: List of parameters of the hemodynamic response function. The list should
%        contain the following:
%        A1 amplitude weight of 1st gamma distribution (default: 0.5)
%        A2 amplitude weight of 2nd gamma distribution (default: 6)
%        A3 amplitude weight of 3rd gamma distribution (default: 1)
%        alpha1 Time to peak of 1st gamma distribution (default: 1.5)
%        alpha2 Time to peak of 2nd gamma distribution (default: 7)
%        alpha3 Time to peak of 3rd gamma distribution (default: 16)
%        beta1 Dispersion of 1st gamma distribution    (default: 0.8)
%        beta2 Dispersion of 2nd gamma distribution    (default: 1)
%        beta3 Dispersion of 3rd gamma distribution    (default: 1)
%        Defaults from Shan et al "Modeling of the hemodynamic responses in block design fMRI
%            studies" J. of Cerebral Blood Flow & Metabolism (2014) 34,
%            316–324.
%
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
    A1 = 0.5;
    A2 = 6;
    A3 = 1;
    alpha1 = 1.5;
    alpha2 = 7;
    alpha3 = 16;    
    beta1 = 0.8;
    beta2 = 1.0;
    beta3 = 1.0;
else
    if ~isstruct(param)
        A1     = param(1);
        A2     = param(2);
        A3     = param(3);
        alpha1 = param(4);
        alpha2 = param(5);
        alpha3 = param(6); 
        beta1  = param(7);
        beta2  = param(8);
        beta3  = param(9);        
    else
        A1     = param.A1;
        A2     = param.A2;
        A3     = param.A3;
        alpha1 = param.alpha1;
        alpha2 = param.alpha2;
        alpha3 = param.alpha3; 
        beta1  = param.beta1;
        beta2  = param.beta2;
        beta3  = param.beta3;          
    end
end

A1_map_flag     = false;
A2_map_flag     = false;
A3_map_flag     = false;
alpha1_map_flag = false;
alpha2_map_flag = false;
alpha3_map_flag = false;
beta1_map_flag  = false;
beta2_map_flag  = false;
beta3_map_flag  = false;
if isfield(A1,'map')
    [~,Y_A1] = STANCE_load_volume(A1.map);    
    A1 = A1.max*double(Y_A1)/255.0;
    A1_map_flag = true;
end
if isfield(A2,'map')
    [~,Y_A2] = STANCE_load_volume(A2.map);    
    A2 = A2.max*double(Y_A2)/255.0;
    A2_map_flag = true;
end
if isfield(A3,'map')
    [~,Y_A3] = STANCE_load_volume(A3.map);    
    A3 = A3.max*double(Y_A3)/255.0;
    A3_map_flag = true;
end
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
if isfield(alpha3,'map')
    [~,Y_alpha3] = STANCE_load_volume(alpha3.map);    
    alpha3 = alpha3.max*double(Y_alpha3)/255.0;
    alpha3_map_flag = true;
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
if isfield(beta3,'map')
    [~,Y_beta3] = STANCE_load_volume(beta3.map);
    beta3 = beta3.max*double(Y_beta3)/255.0;
    beta3_map_flag = true;
end

%% return the triple Gamma function HRF 

if (A1_map_flag) || (A2_map_flag) || (A3_map_flag) || (alpha1_map_flag) || (alpha2_map_flag) || (alpha3_map_flag) || (beta1_map_flag) || (beta2_map_flag) || (beta3_map_flag)
% assume 3D parameter maps
    for i = 1:length(t)
       out(:,:,:,i) = A1.*((t(i).^(alpha1-1)).*(beta1.^alpha1).*exp(-beta1*t(i))./gamma(alpha1)) ...
                   + A2.*((t(i).^(alpha2-1)).*(beta2.^alpha2).*exp(-beta2*t(i))./gamma(alpha2)) ...
                   + A3.*((t(i).^(alpha3-1)).*(beta3.^alpha3).*exp(-beta3*t(i))./gamma(alpha3)); %#ok<AGROW>
    end
else 
    out = A1.*((t.^(alpha1-1)).*(beta1.^alpha1).*exp(-beta1*t)./gamma(alpha1)) ...
        + A2.*((t.^(alpha2-1)).*(beta2.^alpha2).*exp(-beta2*t)./gamma(alpha2)) ...
        + A3.*((t.^(alpha3-1)).*(beta3.^alpha3).*exp(-beta3*t)./gamma(alpha3));    
end

end