function [V_MNI_reg, Y_MNI_reg] = STANCE_register_MNI(filename,M,flags,This_sss)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% Generates a filepath based upon the currently selected study, subject
% and session or This_sss.
%
% Jason E. Hill
% STANCE_register_MNI.m      updated     2 APR 2017

if nargin < 3
    flags = struct('mask',1,'mean',0,'interp',5,'which',1,'wrap',[1 1 0]');
end

if ~exist('STANCE.mat','file')
    if ~exist('..\STANCE.mat','file')
       load('..\..\STANCE.mat'); 
    else
       load('..\STANCE.mat');
    end
else
    load('STANCE.mat');
end

if nargin < 4
    This_sss = Now_sss;
end

[V_MNI,Y_MNI] = STANCE_load_volume(filenameMNI);

% make MNI volume registered to BrainWEB dimensions 
V_MNI_local = V_MNI;
Y_MNI_local = Y_MNI;
Y_max = max(Y_MNI(:));
V_MNI_local.fname = [STANCE_genpath(This_sss,2),'/',filenameMNI(end-18:end-3)]; %#ok<*COLND>
spm_write_vol(V_MNI_local, Y_MNI_local);

spm_get_space(V_MNI_local.fname, M*V_MNI_local.mat);
V_MNI_local = spm_vol(V_MNI_local.fname);

% for re-slicing construct flags first:
% then re-slice, source file first
files = {filename;V_MNI_local.fname};
spm_reslice(files, flags);

delete(V_MNI_local.fname);

fn_MNI_reg = [STANCE_genpath(This_sss,2),'/r',filenameMNI(end-18:end-3)];

[V_MNI_reg,Y_MNI_reg] = STANCE_load_volume(fn_MNI_reg);

% remove negative values (due to 5th degree b-splines)
Y_MNI_reg(Y_MNI_reg < 0) = 0;
Y_MNI_reg(Y_MNI_reg > Y_max) = Y_max;
spm_write_vol(V_MNI_reg, Y_MNI_reg);


