function out = cardiacRF(t, varargin)
% Specifies a double-gamma variate cardiac response function for the 
% given heart rate (HR) time vector and parameters.

%   author:        Dr. Jason E. Hill (post-doc fellow with CNT at TTNI)
%   updated:       24 MAR 2016
%
% Arguments
% t: Time vector in seconds.
% param: List of parameters of the hemodynamic response function. The list should
%        contain the following:
%        a1 Delay of response relative to onset   (default:2.7)
%        b1 Dispersion of response                (default:1.6)
%        b2 Dispersion of undershoot              (default:9)
%        c1 Scale of positive response            (default:0.6)
%        c2 Scale of undershoot                   (default:16)
%        d2 Offset of undershoot                  (default:12)
%         K constant of derivative of CRF_0       (default:0.0)
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
if isempty(param) || (length(param) ~= 5)
    if verbose == true
        warning('Default parameters for HRF are used');
    end
    a1 = 2.7;
    b1 = 1.6;
    c1 = 0.6;
    b2 = 9.0;
    c2 = 16;
    d2 = 12;
    K  = 0;
else
    a1 = param(1);
    b1 = param(2);
    b2 = param(3);
    c1 = param(4);
    c2 = param(5); 
    d2 = param(6);
    K  = param(7);
end

a1_map_flag = false;
b1_map_flag = false;
b2_map_flag = false;
c1_map_flag = false;
c2_map_flag = false;
d2_map_flag = false;
K_map_flag  = false;
if isfield(a1,'map')
    [~,Y_a1] = STANCE_load_volume(a1.map);    
    a1 = a1.max*double(Y_a1)/255.0;
    a1_map_flag = true;
end
if isfield(b1,'map')
    [~,Y_b1] = STANCE_load_volume(b1.map);    
    b1 = b1.max*double(Y_b1)/255.0;
    b1_map_flag = true;
end
if isfield(b2,'map')
    [~,Y_b2] = STANCE_load_volume(b2.map);    
    b2 = b2.max*double(Y_b2)/255.0;
    b2_map_flag = true;
end
if isfield(c1,'map')
    [~,Y_c1] = STANCE_load_volume(c1.map);    
    c1 = c1.max*double(Y_c1)/255.0;
    c1_map_flag = true;
end
if isfield(c2,'map')
    [~,Y_c2] = STANCE_load_volume(c2.map);    
    c2 = c2.max*double(Y_c2)/255.0;
    c2_map_flag = true;
end
if isfield(d2,'map')
    [~,Y_d2] = STANCE_load_volume(d2.map);    
    d2 = d2.max*double(Y_d2)/255.0;
    d2_map_flag = true;
end
if isfield(K,'map')
    [~,Y_K] = STANCE_load_volume(K.map);    
    K = K.max*double(Y_K)/255.0;
    K_map_flag = true;
end

%% return the CRF
if (a1_map_flag) || (b1_map_flag) || (b2_map_flag) || (c1_map_flag) || (c2_map_flag) || (d2_map_flag)  || (K_map_flag)
% assume 3D parameter maps
    for i = 1:length(t)
       DtCRF = c1.*(t.^a1).*exp(-t./b1).*(a1./t - 1./b1) - ((t-d2)./b2).*(c2./sqrt(2*pi.*b2)).*exp(-0.5.*(t-d2).^2./b2);
       out(:,:,:,i) = c1.*(t(i).^a1).*exp(-t(i)./b1) - (c2./sqrt(2*pi.*b2)).*exp(-0.5.*(t(i)-d2).^2./b2) + K.*DtCRF; %#ok<AGROW>
    end
else
    DtCRF = c1.*(t.^a1).*exp(-t./b1).*(a1./t - 1./b1) - ((t-d2)./b2).*(c2./sqrt(2*pi.*b2)).*exp(-0.5.*(t-d2).^2./b2);
    out = c1.*(t.^a1).*exp(-t./b1) - (c2./sqrt(2*pi.*b2)).*exp(-0.5.*(t-d2).^2./b2) + K.*DtCRF;    
end
% Reference:
% Dietmar Cordesa, Rajesh R. Nandy, Scott Schafer, and Tor D. Wager,
% "Characterization and Reduction of Cardiac- and Respiratory-Induced Noise as a Function of the Sampling Rate (TR) in fMRI"
% Neuroimage. 2014 April 1; 89:314-330.

end