function M_array = STANCE_find_affine_array
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% prepares affine tranformation matrices M for use with STANCE
% s is the subject number

if ~exist('STANCE.mat','file')
    if ~exist('..\STANCE.mat','file')
       load('..\..\STANCE.mat'); 
    else
       load('..\STANCE.mat');
    end
else
    load('STANCE.mat');
end

M_array  = zeros(4,4,20);        
        
for s = 1:20
    dispString = sprintf('o Computing transform for subject %d.',s);
    disp(dispString)
    
    if s < 4
       subject  = [subjectBase,'0', num2str(subjects(s))];
    else
       subject  = [subjectBase, num2str(subjects(s))];       
    end
    fileNameNII  = [subject,'_t1w_p0.nii.gz'];
    
    [V_s,~] = STANCE_load_volume(fileNameNII);
    
    [V_MNI,~] = STANCE_load_volume(filenameMNI);

    x = spm_coreg(V_MNI,V_s);
    M = spm_matrix(x);
    
    M_array(:,:,s) = M;
    
end

%M_std = std(M_array,[],3);
currentDir = pwd;
cd(STANCEroot)
save('STANCE.mat','SPMpath','subject_labels','Nbr_subject_labels','subjectsBase','subjectsT1postfix','subjectsFuzzyPostfix','filenameMNI','filenameGM','Nbr_sss','Now_sss','STANCEroot','male_labels','female_labels','M_array','MNIgmVolume');
cd(currentDir)    
