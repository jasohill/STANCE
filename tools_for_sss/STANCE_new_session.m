function [Nbr_sss,Now_sss] = STANCE_new_session(study_number,subject_number,new_session_number,OKflag)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE)
%
% Creates file paths for new session and updates global variables.
%
% Jason E. Hill
% STANCE_new_session.m      updated     2 APR 2017

if ~exist('STANCE.mat','file')
    load('..\STANCE.mat');
else
    load('STANCE.mat');
end
oldDir = pwd;

if nargin < 4
    OKflag = false;
end

if nargin < 1
    study_number = Now_sss(1); %#ok<*NODEF>
end

no3arginFlag = 0;
if nargin < 2
    if length(study_number) == 3
        This_sss = study_number;
        study_number = This_sss(1); 
        subject_number = This_sss(2);    
        new_session_number = This_sss(3);
        no3arginFlag = 1;     
    else
        subject_number = Now_sss(2);
    end
end

if nargin < 3 && ~no3arginFlag
    Nbr_sessions = STANCE_how_many_sessions(study_number,subject_number);    
    new_session_number = Nbr_sessions + 1;
    cd(oldDir);
elseif isempty(new_session_number) || new_session_number == 0 
    if new_session_number == 0
        Nbr_sessions = 0;
    else
        Nbr_sessions = STANCE_how_many_sessions(study_number,subject_number);
    end
    new_session_number = Nbr_sessions + 1;
    cd(oldDir);    
end

Temp_sss(1) = study_number;
Temp_sss(2) = subject_number;
Temp_sss(3) = new_session_number;

stopFlag = false;
filepath = STANCE_genpath(Temp_sss);

if logical(isdir(filepath)) && ~stopFlag;
    if length(dir(filepath)) > 2
        % Construct a questdlg with two options
        query = ['This folder is not empty! Remove files in ' filepath '?'];
        choice.default = 'OK'; %#ok<STRNU>
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
    end
else
    mkdir(filepath);
end

if stopFlag
    % do not update anything
else
    Now_sss = Temp_sss;
    if sum(Now_sss)>3
        Nbr_sss(3) = Nbr_sss(3) + 1;   
    end
    currentDir = pwd;
    cd(STANCEroot)
    save('STANCE.mat','SPMpath','subject_labels','Nbr_subject_labels','subjectsBase','subjectsT1postfix','subjectsFuzzyPostfix','filenameMNI','filenameGM','Nbr_sss','Now_sss','STANCEroot','male_labels','female_labels','M_array','MNIgmVolume');
    cd(currentDir)    
end
