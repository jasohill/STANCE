function Nbr_subjects = STANCE_how_many_subjects(study_number)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% Returns the number of subjects in an study
%
% Jason E. Hill
% STANCE_how_many_subjects.m      updated     2 OCT 2016

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
        list = ls('subject*');
        Nbr_subjects = size(list,1);
        cd(oldDir);
    else
        Nbr_subjects = 0;
    end
end

