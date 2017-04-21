function [Nbr_sss,Now_sss] = STANCE_new_study(new_study_number,OKflag)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%
% Creates file paths for new study and updates global variables.
%
% Jason E. Hill
% STANCE_new_study.m      updated     2 APR 2017

if ~exist('STANCE.mat','file')
    load('..\STANCE.mat');
else
    load('STANCE.mat');
end

if nargin < 2
    OKflag = false;
end

if nargin < 1
    [Nbr_studies,~] = STANCE_how_many_studies;
    This_sss(1) = Nbr_studies + 1;
    This_sss(2) = 1;
    This_sss(3) = 1;    
    new_study_number = This_sss(1);
end

This_sss(1) = new_study_number;
This_sss(2) = 1;
This_sss(3) = 1;

Temp_sss(1) = This_sss(1);
Temp_sss(2) = This_sss(2);
Temp_sss(3) = This_sss(3);

stopFlag = false;
for  level = 1:3
    
    filepath = STANCE_genpath(Temp_sss,level);
    
    if isdir(filepath) && ~stopFlag; 
        % Construct a questdlg with two options
        query = ['This folder is not empty! Remove files in ' filepath '?'];
        choice.default = 'OK';
        if OKflag
            choice = 'OK';
        else
            choice = questdlg(query, ...
                         'YesNo Menu', ...
                         'OK','No','No');
        end
        % Handle response
        switch choice
            case 'OK'
                disp(['Erasing files in ' filepath '.'])
                cmd_rmdir(filepath);
                mkdir(filepath);
            case 'No.'
                disp(['Cannot create new study in' filepath '. One already exists there.'])
                stopFlag = true;
        end
    else        
        mkdir(filepath);
    end
end

if stopFlag
    % do not update anything
else
    Now_sss = Temp_sss;
    Nbr_sss(1) = Nbr_sss(1) + 1;    
    currentDir = pwd;
    cd(STANCEroot)
    save('STANCE.mat','SPMpath','subject_labels','Nbr_subject_labels','subjectsBase','subjectsT1postfix','subjectsFuzzyPostfix','filenameMNI','filenameGM','Nbr_sss','Now_sss','STANCEroot','male_labels','female_labels','M_array','MNIgmVolume');
    cd(currentDir)    
end

