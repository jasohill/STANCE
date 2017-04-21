function session_list = STANCE_list_sessions(study_number,subject_number)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% Returns the number of subjects in an study
%
% Jason E. Hill
% STANCE_list_sessions.m      updated     2 OCT 2016

if ~exist('STANCE.mat','file')
    load('../STANCE.mat');
else
    load('STANCE.mat');
end

    if nargin < 1
        study_number = Now_sss(1);
    end
    if nargin < 2
        subject_number = Now_sss(2);
    end
    subjectFilepath = STANCE_genpath([study_number,subject_number,1],2);

    if isdir(subjectFilepath)
        oldDir = cd(subjectFilepath);    
        session_list = ls('session*');
        cd(oldDir);
    else
        session_list = [];
    end
end