function subject_list = STANCE_list_subjects(study_number)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% Returns the number of subjects in an study
%
% Jason E. Hill
% STANCE_list_subjects.m      updated     2 OCT 2016

if ~exist('STANCE.mat','file')
    load('../STANCE.mat');
else
    load('STANCE.mat');
end

    if nargin < 1
        study_number = Now_sss(1);
    end

    studyFilepath = STANCE_genpath([study_number,1,1],1);

    if isdir(studyFilepath)
        oldDir = cd(studyFilepath);    
        subject_list = ls('subject*');
        cd(oldDir);
    else
        subject_list = [];
    end
end