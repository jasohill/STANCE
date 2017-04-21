function filepath = STANCE_genpath(This_sss,level)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% Generates a filepath based upon the currently selected study, subject
% and session
%
% Jason E. Hill
% STANCE_genpath.m      updated     18 FEB 2016

if nargin < 2
    level = 3;
end

if nargin < 1
    This_sss = [];
end

if ~exist('STANCE.mat','file')
    load('../STANCE.mat');
else
    load('STANCE.mat');
end

if isempty(This_sss)
    This_sss = Now_sss;
end
study_number   = This_sss(1);
subject_number = This_sss(2);
session_number = This_sss(3);

study   = sprintf('%s%.3d','study'  ,study_number);  
subject = sprintf('%s%.4d','subject',subject_number);  
session = sprintf('%s%.3d','session',session_number);  

if level > 2
    filepath  = [STANCEroot,'/','fMRI/',study,'/',subject,'/',session]; 
elseif level == 2
    filepath  = [STANCEroot,'/','fMRI/',study,'/',subject];     
else
    filepath  = [STANCEroot,'/','fMRI/',study];     
end

end

