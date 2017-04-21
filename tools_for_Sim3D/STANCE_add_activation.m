function [V_act,Y_act] = STANCE_add_activation(filename_baseline,activation_map,TE,amplitude)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%
% Adds activation to T2* baseline map to generate BOLD signal.
%
% Jason E. Hill, Ph.D.
% STANCE_add_activation.m      updated     27 MAR 2017

if nargin <4
    amplitude = 0.02;  % default is 2% of signal
end

if nargin < 3
    TE = 30; % [ms]
end

[V_baseline,Y_baseline] = STANCE_load_volume(filename_baseline);

% add activation to T2* baseline
Delta_S = amplitude.*activation_map;  
Y_act = Y_baseline.*(1 + Delta_S.*(Y_baseline./TE).*(1+((Y_baseline./TE)-0.5).*Delta_S));
V_act = V_baseline;
[fileBase,fileExtension] = get_file_parts(filename_baseline);
V_act.fname = [fileBase,'_BOLD',fileExtension];
V_act = spm_write_vol(V_act,Y_act);

end

