function [V_act_reg,Y_act_reg] = STANCE_register_activation(filename,task,flags,This_sss,ADO)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% Registers task.map to filename's volume and saves in the subject of interest 
% folder as a NIFTI file named ‘r’+task.name.
%
% Jason E. Hill
% STANCE_register_activation.m      updated     3 APR 2017

if nargin < 5
    ADO = 1; % assume anatomical space
end
if nargin < 3
    flags = struct('interp',5,'mask',1,'mean',0,'which',1,'wrap',[1 1 0]');
end
if isempty(flags)
    flags = struct('interp',5,'mask',1,'mean',0,'which',1,'wrap',[1 1 0]');
end    

if ~exist('STANCE.mat','file')
    if ~exist('../STANCE.mat','file')
       load('../../STANCE.mat'); 
    else
       load('../STANCE.mat');
    end
else
    if exist('STANCE.mat','file')
       load('STANCE.mat');
    end
end
if exist(SPMpath,'dir')
    addpath(SPMpath);    
end
if nargin < 4 || isempty(This_sss)
    This_sss = Now_sss;
end

if ADO == 1
    activationMap = task.map;
else
    activationMap = task.activation(ADO-1).map;
end    

[V_MNI,~] = STANCE_load_volume(filenameMNI);

V_act = V_MNI;
V_act.fname = [STANCEroot,'/activations/',task.name,'.nii'];
V_act = spm_create_vol(V_act);
V_act = spm_write_vol(V_act,activationMap);

files = {filename;V_act.fname};
spm_reslice(files, flags);
fn_act_reg = [STANCEroot,'/activations/r',task.name,'.nii'];
[V_act_reg,Y_act_reg] = STANCE_load_volume(fn_act_reg);
delete(fn_act_reg);
fn_act_reg = [STANCE_genpath(This_sss,2),'/r',task.name,'.nii'];
V_act_reg.fname = fn_act_reg;

% remove negative values (due to 5th degree b-splines)
Y_act_reg(Y_act_reg < 0) = 0;
Y_act_reg(Y_act_reg > 1.0) = 1.0;


V_act_reg = spm_write_vol(V_act_reg,Y_act_reg);

[V_act_reg,Y_act_reg] = STANCE_load_volume(fn_act_reg);

delete(V_act.fname);

end
