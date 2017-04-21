function M_array = prep_subject_files
% Prepares subject files for use with STANCE by unzipping them and
% registering the T1-w MNI152 normalized brain to each of them.
%
% Jason E. Hill
% prep_subject_files.m      updated     2 APR 2017

subjects = [4,5,6,18,20,38,41,42,43,44,...
            45,46,47,48,49,50,51,52,53,54];

M_array  = zeros(4,4,12);        
        
for s = 1:20 % s is the subject number
    if s < 4
       subject  = ['subject0', num2str(subjects(s))];
    else
       subject  = ['subject', num2str(subjects(s))];       
    end        
    fileNameNIIt1  = [subject,'_t1w_p0.nii.gz'];     
    fileNameNIIt1Out  = ['subjects\',subject,'\',subject,'_t1w_p0.nii'];     
    fileNameNIIfuzzy = [subject,'_fuzzy.nii.gz']; 
    fileNameNIIfuzzyOut  = ['subjects\',subject,'\',subject,'_fuzzy.nii']; 
   
    [V_sXX,Y_sXX] = STANCE_load_volume(fileNameNIIt1);

    fn_MNI = 'MNI\MNI152_T1_1mm_181x217x181.nii.gz';
    [V_MNI,~] = STANCE_load_volume(fn_MNI);

    x = spm_coreg(V_MNI,V_sXX);
    M = spm_matrix(x);

    V_sXX.private.mat = M;
    V_sXX.private.mat0_intent = 'UNKNOWN'; % used to store variablity models (temp)
    V_sXX.private.mat0 = M;
    V_sXX.descrip = 'BrainWEB 20 Norm';
    V_sXX.fname = fileNameNIIt1Out;
    V_sXX.private.mat_intent = 'MNI152';

    spm_write_vol(V_sXX, Y_sXX);

    gzip(V_sXX.fname);
    delete(V_sXX.fname);

    [V_sXX_tissue,Y_sXX_tissue] = STANCE_load_volume(fileNameNIIfuzzy);

    header = V_sXX_tissue(1);
    header.descrip = 'BrainWEB 20 Norm';
    header.fname = fileNameNIIfuzzyOut;
    header.private.mat_intent = 'MNI152';
    header.private.mat = M;
    header.private.mat0_intent = 'UNKNOWN'; % used to store variablity models (temp)
    header.private.mat0 = M;
    header.private.timing.tspace = 0;

    for vol = 1:12
       header.n(1) = vol;    
       spm_write_vol(header,Y_sXX_tissue(:, :, :, vol));
    end

    gzip(header.fname);
    delete(header.fname);

    M_array(:,:,s) = M;
end

end

