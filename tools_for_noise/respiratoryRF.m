function out = respiratoryRF(t, varargin)
% Specifies a double-gamma variate respiratory response function for the 
% given respiratory volume time-series (RVT) time vector and parameters.

%   author:         Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
%   updated:        16 NOV 2016
%
% Arguments
% t: Time vector in seconds.
% param: List of parameters of the respiratory response function. The list should
%        contain the following:
%        a1 Delay of response relative to onset   (default:2.1)
%        a2 Delay of undershoot relative to onset (default:3.54)
%        b1 Dispersion of response                (default:4.25)
%        b2 Dispersion of undershoot              (default:0.9)
%        c1 Scale of undershoot                   (default:0.6)
%        c2 Scale of undershoot                   (default:0.0023)
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
    
    a1 = 2.1;
    b1 = 1.6;
    c1 = 0.6;    
    a2 = 3.54;
    b2 = 4.25;
    c2 = 0.0023;
else
    a1 = param(1);
    a2 = param(2);
    b1 = param(3);
    b2 = param(4);
    c1 = param(5); 
    c2 = param(6);     
end

a1_map_flag = false;
a2_map_flag = false;
b1_map_flag  = false;
b2_map_flag  = false;
c1_map_flag  = false;
c2_map_flag  = false;
if isfield(a1,'map')
    [~,Y_a1] = STANCE_load_volume(a1.map);    
    a1 = a1.max*double(Y_a1)/255.0;
    a1_map_flag = true;
end
if isfield(a2,'map')
    [~,Y_a2] = STANCE_load_volume(a2.map);    
    a2 = a2.max*double(Y_a2)/255.0;
    a2_map_flag = true;
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

%% return the double Gamma function RRF

if (a1_map_flag) || (a2_map_flag) || (b1_map_flag) || (b2_map_flag) || (c1_map_flag) || (c2_map_flag)
    % assume 3D parameter maps
    for i = 1:length(t)
       out(:,:,:,i) = c1.*(t(i).^a1).*exp(-t(i)./b1) - c2.*(t(i).^a2).*exp(-t(i)./b2); %#ok<AGROW>
    end
else 
    out = c1*(t.^a1).*exp(-t./b1) - c2.*(t.^a2).*exp(-t./b2);
end

end