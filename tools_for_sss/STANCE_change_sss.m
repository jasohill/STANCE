function Now_sss = STANCE_change_sss(session_number,subject_number,study_number)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% Switch to another study, subject, session.
%
% Jason E. Hill
% STANCE_change_sss.m      updated     18 FEB 2016


if ~exist('STANCE.mat','file')
    load('..\STANCE.mat');
else
    load('STANCE.mat');
end

if nargin < 3
    study_number = Now_sss(1); %#ok<NODEF>
end

if nargin < 2
    subject_number = Now_sss(2);
end

if nargin < 1
    session_number = Now_sss(3) + 1;
end

Temp_sss = [study_number, subject_number, session_number];

filepath = STANCE_genpath(Temp_sss);

if ~isdir(filepath)
   warning('Cannot comply. File path does not exist for this selection!') 
else
    Now_sss = Temp_sss;
    currentDir = pwd;
    cd(STANCEroot)
    save('STANCE.mat','SPMpath','subject_labels','Nbr_subject_labels','subjectsBase','subjectsT1postfix','subjectsFuzzyPostfix','filenameMNI','filenameGM','Nbr_sss','Now_sss','STANCEroot','male_labels','female_labels','M_array','MNIgmVolume');
    cd(currentDir)
end
    
end

