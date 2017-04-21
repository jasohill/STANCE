function [V_out, Y_out] = STANCE_conform(P_in, P_ref, file_out, description, gzipFlag)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% Conforms V_in to the dimensions with same origin and orientation of V_ref. 
% Negative diagonal elements of mat encodes orientation.
%
% Jason E. Hill
% STANCE_conform.m      updated     14 SEPT 2016

if nargin < 5
    gzipFlag = false; % do NOT save output file as a *.gz archive file
end

if nargin < 4
    description = [];
end

if ischar(P_in)
    V_in = STANCE_load_header(P_in);
elseif isstruct(P_in)
    V_in = P_in;    
end

if ischar(P_ref)
    V_ref = STANCE_load_header(P_ref);
elseif isstruct(P_ref)
    V_ref = P_ref;    
end

ref_dimensions = V_ref(1).dim;
ref_origin = round(abs(V_ref(1).mat(1:3,4)));
ref_mat = V_ref(1).mat;
ref_orient = sign([ref_mat(1,1), ref_mat(2,2), ref_mat(3,3)]);  

in_size = size(V_in); %#ok<NASGU>
in_fn = V_in(1).fname;
in_dimensions = V_in(1).dim;
in_mat = V_in(1).mat; 
in_origin = round(abs(in_mat(1:3,4)));
in_orient = sign([in_mat(1,1), in_mat(2,2), in_mat(3,3)]);  

Y_in = STANCE_load_data(in_fn);

origin_diff = (ref_origin - in_origin);
orient_flip = (ref_orient.*in_orient);

out_mat = in_mat;
for d = 1:3
    out_mat(d,d) = sign(in_mat(d,d))*abs(ref_mat(d,d));
end
%% transform scale
if (in_mat(1,1) ~= ref_mat(1,1))||(in_mat(2,2) ~= ref_mat(1,1))||(in_mat(3,3) ~= ref_mat(3,3))
    for d = 1:3
       scale_dimensions(d) = abs(ceil(in_dimensions(d)*(in_mat(d,d)/ref_mat(d,d)))); %#ok<AGROW>
    end
    reslice(in_fn,file_out,scale_dimensions,out_mat,7);
    in_dimensions = scale_dimensions;
    
    Y_max = max(Y_in(:));

    Y_in = STANCE_load_data(file_out);
    
    % correct limits
    Y_in(Y_in<0) = 0.0;
    Y_in(Y_in>Y_max) = Y_max;
else
    Y_in = STANCE_load_data(in_fn);
end

dimens_diff = (ref_dimensions - in_dimensions);

%% transform origin
for d = 1:3
    out_mat(d,d) = ref_mat(d,d);
    out_mat(d,4) = sign(ref_mat(d,4))*(abs(in_mat(d,4)) + origin_diff(d));
end

padding  = double(uint8(dimens_diff));
postpad = double(uint8(in_orient.*padding));
prepad = double(uint8(-in_orient.*padding));

Y_pad = padarray(padarray(Y_in,prepad,0,'pre'),postpad,0,'post');

chopping = double(uint8(-dimens_diff));
postchop = uint8(in_orient.*chopping);
prechop = uint8(-in_orient.*chopping);

Y_out = zeros(ref_dimensions); %#ok<PREALL>
Y_out = Y_pad(1+prechop(1):max([in_dimensions(1),ref_dimensions(1)])-postchop(1),1+prechop(2):max([in_dimensions(2),ref_dimensions(2)])-postchop(2),1+prechop(3):max([in_dimensions(3),ref_dimensions(3)])-postchop(3));

for d = 1:3
   if orient_flip(d) == -1
     Y_out = flip(Y_out,d);  
   end
end

Y_out = circshift(Y_out,origin_diff);

V_out = V_in;
V_out(:).fname = file_out;
V_out(:).dim = ref_dimensions;
V_out(:).mat = out_mat;
if ~isempty(description)
    V_out(:).descrip = description;
end

V_out = spm_create_vol(V_out);
V_out = spm_write_vol(V_out,Y_out);

if gzipFlag
    ungzipFilename = file_out;
    file_out = gzip(file_out);
    delete(ungzipFilename);
    if iscellstr(file_out)
        file_out = char(file_out{1,1});
    end
end
V_out(:).fname = char(file_out);

end

