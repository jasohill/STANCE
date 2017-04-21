function map = STANCE_load_map(P_in, P_ref, specification, gzipFlag)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% load P_in as a map to use as an activation region 
% with the same dimensions and origin as P_ref
%   
%   gzipFlag indicated whether to keep only archived gzipped files (lower 
%   storage cost) or allowing for ungzipped files to remain (faster run times).
%
% author: Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
% STANCE_load_map.m      updated     24 MAR 2016

if ~exist('STANCE.mat','file')
    if ~exist('..\STANCE.mat','file')
       load('..\..\STANCE.mat'); 
    else
       load('..\STANCE.mat');
    end
else
    load('STANCE.mat');
end

if nargin < 4
    gzipFlag  = true;
end

if nargin < 3
    specification = [];
end

if ischar(P_in)
    fn_in = P_in;
elseif isstruct(P_in)
    V_in = P_in;    
    fn_in = V_in.fname;
end
[V_in,Y_in] = STANCE_load_volume(fn_in,[],gzipFlag);

if ischar(P_ref)
    V_ref = STANCE_load_header(P_ref);
elseif isstruct(P_ref)
    V_ref = P_ref;    
end

if ndims(Y_in) > 3
    if isempty(specification)
        specification = 1;
    end
    msg = ['o Loading map from component data in ',fn_in,'.'];
    disp(msg)
    if size(Y_in,4) > 1      
        if round(abs(specification)) == 0
            Y_act = squeeze(Y_in(:,:,:,1));           
        else
            specification = round(abs(specification));  
            Y_act = squeeze(Y_in(:,:,:,specification));
        end
    else
        Y_act = squeeze(Y_in(Y_in>specification));
    end
    V_act = V_in(1);
else
    if isempty(specification)
        threshold = 0; %#ok<NASGU>
        saturation = max(Y_in(:)); %#ok<NASGU>
    end    
    if length(specification) == 1
        threshold = specification;
        saturation = max(Y_in(:));        % 6 sigma standard        
    else
        threshold = specification(1);  % lower threshold for activation
        saturation = specification(2); % upper threshold for saturation       
    end
    msg = ['o Loading map from data in ',fn_in,'.'];
    disp(msg)

    Y_act = Y_in.*(Y_in > threshold);
    Y_act(Y_act > saturation) = saturation;
    Y_act = Y_act./saturation;
    V_act = V_in;
end
file_act_in  = [STANCEroot,'/Activations/temp_act.nii'];
file_act_out = [STANCEroot,'/Activations/temp_act_1mm.nii'];
V_act.fname = file_act_in;
V_act = spm_create_vol(V_act);
V_act = spm_write_vol(V_act,Y_act);

[V_out,Y_out] = STANCE_conform(V_act,V_ref,file_act_out,'Resliced to MNI 1mm^3'); %#ok<ASGLU>

% normalize range to 0..1
minActOut = min(Y_out(:));
maxActOut = max(Y_out(:));
map = (Y_out + minActOut)./(maxActOut-minActOut);

if ~strcmp(fn_in,V_in(1).fname) && gzipFlag
    delete(V_in(1).fname)
end
delete(file_act_in);
delete(file_act_out);

end