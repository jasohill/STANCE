function STANCE_initialize_STANCE
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%
% Initializes STANCE.mat and other global resources, such as many of the file 
% locations, such as the provided 20 subject brains and their file names, 
% upon the first execution of STANCE
%
% Jason E. Hill
% STANCE_initialize_STANCE.m      updated     17 FEB 2016
%clear all, close all;  %#ok<CLALL>

addpath(genpath(pwd));

display('Welcome to STANCE, starting initialization...')

% if exist(spm('Dir'),'dir')
%     display('o SPM installation found.')
%     SPMpath = spm('Dir');
% else
    hSPM8 = msgbox('Please select the SPM8 directory');
    uiwait(hSPM8);
    currPath = fileparts(mfilename('fullpath'));
    SPMpath = uigetdir(currPath, 'Add SPM filepath');
    addpath(genpath(pwd));
    if SPMpath
        if strcmp(SPMpath(:,end-3:end),'spm8')
            addpath(SPMpath);
        else
            hSPM8=errordlg('SPM path is wrong!','ERROR');
            ha=get(hSPM8,'children');
            hu=findall(allchild(hSPM8),'style','pushbutton');
            set(hu,'string','OK');
            ht=findall(ha,'type','text');
            set(ht,'fontsize',12);
            movegui(hSPM8, 'center');
        end
    else
        return
    end
% end

%% Turn off ...
% ... OpenGl warnings
warning('off','MATLAB:opengl:StartupBlacklistedNoSetting');
warning('off', 'MATLAB:hg:AutoSoftwareOpenGL');
% ... finite warning
warning('off', 'MATLAB:FINITE:obsoleteFunction');
% ... NIFTI class warnings when loading SPM mat files
warning('off', 'MATLAB:unknownElementsNowStruc');
warning('off', 'MATLAB:dispatcher:ShadowedMEXExtension');
warning('off', 'MATLAB:pfileOlderThanMfile');
% ... removing files from path
warning('off', 'MATLAB:RMDIR:RemovedFromPath');


%% STANCE global variables
% determine the root STANCE directory
currentDir = pwd;
if strcmp(currentDir(end-5:end),'STANCE')
    STANCEroot = pwd;    
elseif strcmp(currentDir(end-2:end),'GUI')
    % GUI instance of initialization
    cd ../
    STANCEroot = pwd;
    cd(currentDir)
else
    % call originates elsewhere must locate the STANCE folder 
    hSTANCE = msgbox('Please the select STANCE directory');
    uiwait(hSTANCE);
    currPath = fileparts(mfilename('fullpath'));
    STANCEroot = uigetdir(currPath, 'Add SPM filepath');
    addpath(genpath(pwd));
    if STANCEroot
        if strcmp(STANCEroot(:,end-5:end),'STANCE')
            addpath(STANCEroot);
        else
            hSTANCE=errordlg('STANCE path is wrong!','ERROR');
            ha=get(hSTANCE,'children');
            hu=findall(allchild(hSTANCE),'style','pushbutton');
            set(hu,'string','OK');
            ht=findall(ha,'type','text');
            set(ht,'fontsize',12);
            movegui(hSTANCE, 'center');
        end
    else
        return
    end    
end
filenameMNI = [STANCEroot,'/MNI/MNI152T1_1mm.nii'];
filenameGM  = [STANCEroot,'/MNI/MNI152gm_1mm.nii'];

%subject related objects
subject_labels = [4,5,6,18,20,38,41,42,43,44,45,46,47,48,49,50,51,52,53,54];
male_labels    = [1,3,4,6,9,12,16,17,19,20]; %#ok<*NASGU>
female_labels  = [2,5,7,8,10,11,13,14,15,18];    
Nbr_subject_labels = length(subject_labels);
Nbr_sss = [1,1,1];   % number of studies, subjects & sessions
Now_sss = [1,1,1];   % study, subject and session being simulated
subjectsBase = [STANCEroot,'/subjects/subject'];
subjectsT1postfix = '_t1w_p0.nii.gz';
subjectsFuzzyPostfix = '_fuzzy.nii.gz';

display('o Generating 1mm MNI files from SPM resources.')

spmMNI = [spm('Dir'),'/canonical/avg152T1.nii'];
spmGM = [spm('Dir'),'/apriori/grey.nii'];
fn_s04 = 'subjects/subject04/subject04_t1w_p0.nii.gz';
V_s04 = STANCE_load_header(fn_s04);

currentDir = cd(STANCEroot);

[V_MNI,~] = STANCE_conform(spmMNI,V_s04,filenameMNI,'Resliced SPM',true);
[~,Y_GM]  = STANCE_conform(spmGM,V_s04,filenameGM,'Resliced SPM',true);
filenameMNI = [filenameMNI,'.gz'];
filenameGM  = [filenameGM,'.gz'];
  
STANCEflag = false;
if exist(fullfile(cd,'STANCE.mat'),'file')
    load('STANCE.mat')
    warning('Attempting to initialize when STANCE.mat already exists!')
    STANCEflag = true;
end

cd(currentDir);

MNIgmVolume = sum(Y_GM(:));

if ~STANCEflag
    
    M_array  = zeros(4,4,12);
    V_MNIu = STANCE_load_header(V_MNI.fname);

    display('o Calculating the affine transformation matrices for 20 brains, please be patient.')
    for s = 1:20
        dispString = sprintf('o Computing transform for subject %d.',s);
        disp(dispString)
        if s < 4
           subject  = ['subject0', num2str(subject_labels(s))];
        else
           subject  = ['subject', num2str(subject_labels(s))];       
        end        
        fileNameNIIt1  = [subjectsBase,subject(8:9),'/',subject,'_t1w_p0.nii.gz'];      
   
        V_sXX = STANCE_load_header(fileNameNIIt1);

        x = spm_coreg(V_MNIu,V_sXX);
        M = spm_matrix(x);
        M_array(:,:,s) = M;
    end
end

if ~STANCEflag
    display('o Saving initialization to STANCE.mat.')
    currentDir = pwd;
    cd(STANCEroot)
    save('STANCE.mat','SPMpath','subject_labels','Nbr_subject_labels','subjectsBase','subjectsT1postfix','subjectsFuzzyPostfix','filenameMNI','filenameGM','Nbr_sss','Now_sss','STANCEroot','male_labels','female_labels','M_array','MNIgmVolume');
    cd(currentDir)
    
    display('o Creating new study, subject and session folders.')
    STANCE_new_session(1,1,0);
end

display('Done with STANCE initialization.')

end

