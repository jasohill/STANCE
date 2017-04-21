function GM_volume = STANCE_find_GM_volume(task,ADO)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% Returns the relative MNI gray matter volume of activation map
%
% Jason E. Hill
% STANCE_find_GM_volume.m      updated     29 SEPT 2016

if nargin < 2
   ADO = 1;
end
ADO = round(ADO);
if ADO < 1
   ADO = 1;
end
if ADO == 1
   activation_map = task.map;
else
   activation_map = task.activation(ADO-1).map;
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

filenameGM = [STANCEroot,'\MNI\MNI152gm_1mm.nii.gz'];

[V_GM,Y_GM] = STANCE_load_volume(filenameGM);
    
GM_volume = sum(sum(sum(Y_GM.*activation_map)))/MNIgmVolume;
