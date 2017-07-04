function out = tripleLogitHRF(t, varargin)
% Specifies a triple-logit variate hemodynamic response function for the 
% given time vector and parameters.

%   author:        Dr. Jason E. Hill (post-doc fellow with CNT at TTNI)
%   updated:       25 MAR 2017
%
% Arguments
% t: Time vector in seconds.
% param: List of parameters of the hemodynamic response function. The list should
%        contain the following:
%        A1 amplitude weight of 1st inverse logit function (default: 1)
%        A2 amplitude weight of 2nd inverse logit function (default: 1)
%        A3 amplitude weight of 3rd inverse logit function (default: 1)
%        T1 shift center of 1st inverse logit function     (default: 4)
%        T2 shift center of 2nd inverse logit function     (default: 5)
%        T3 shift center of 3rd inverse logit function     (default: 10)
%        D1 slope of 1st inverse logit function            (default: 1)
%        D2 slope of 2nd inverse logit function            (default: 1.5)
%        D3 slope of 3rd inverse logit function            (default: 2)
%        Defaults from Shan et al "Modeling of the hemodynamic responses in block design fMRI
%            studies" J. of Cerebral Blood Flow & Metabolism (2014) 34(2),
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
    A1 = 1;
    A2 = 1;
    A3 = 1;    
    T1 = 4;
    T2 = 5;
    T3 = 10;    
    D1 = 1;
    D2 = 1.5;
    D3 = 2;
else
    if ~isstruct(param)
        A1     = param(1);
        A2     = param(2);
        A3     = param(3);        
        T1     = param(4);
        T2     = param(5);
        T3     = param(6);
        D1     = param(7);
        D2     = param(8);
        D3     = param(9);        
    else
        A1     = param.A1;
        A2     = param.A2;
        A3     = param.A3;        
        T1     = param.T1;
        T2     = param.T2;
        T3     = param.T3;
        D1     = param.D1;
        D2     = param.D2;
        D3     = param.D3;          
    end
end

A1_map_flag     = false;
A2_map_flag     = false;
A3_map_flag     = false;
T1_map_flag     = false;
T2_map_flag     = false;
T3_map_flag     = false;
D1_map_flag     = false;
D2_map_flag     = false;
D3_map_flag     = false;
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
if isfield(T1,'map')
    [~,Y_T1] = STANCE_load_volume(T1.map);    
    T1 = T1.max*double(Y_T1)/255.0;
    T1_map_flag = true;
end
if isfield(T2,'map')
    [~,Y_T2] = STANCE_load_volume(T2.map);    
    T2 = T2.max*double(Y_T2)/255.0;
    T2_map_flag = true;
end
if isfield(T3,'map')
    [~,Y_T3] = STANCE_load_volume(T3.map);    
    T3 = T3.max*double(Y_T3)/255.0;
    T3_map_flag = true;
end
if isfield(D1,'map')
    [~,Y_D1] = STANCE_load_volume(D1.map);    
    D1 = D1.max*double(Y_D1)/255.0;
    D1_map_flag = true;
end
if isfield(D2,'map')
    [~,Y_D2] = STANCE_load_volume(D2.map);    
    D2 = D2.max*double(Y_D2)/255.0;
    D2_map_flag = true;
end
if isfield(D3,'map')
    [~,Y_D3] = STANCE_load_volume(D3.map);    
    D3 = D3.max*double(Y_D3)/255.0;
    D3_map_flag = true;
end

%% return the triple logit function HRF 

if (A1_map_flag) || (A2_map_flag) || (A3_map_flag) || (T1_map_flag) || (T2_map_flag) || (T3_map_flag) || (D1_map_flag) || (D2_map_flag) || (D3_map_flag)
    % assume 3D parameter maps
    for i = 1:length(t)
       out(:,:,:,i) = A1./(1 + exp((t(i)-T1)./D1)) ...
                    + A2./(1 + exp((t(i)-T2)./D2)) ...
                    + A3./(1 + exp((t(i)-T3)./D3)); %#ok<AGROW>
    end
else 
    out = A1./(1 + exp((t-T1)./D1)) ...
        + A2./(1 + exp((t-T2)./D2)) ...
        + A3./(1 + exp((t-T3)./D3));
end

end