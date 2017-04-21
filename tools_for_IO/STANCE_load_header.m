function V = STANCE_load_header(fileName,dataDir,gzipFlag)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%STANCE_load_header loads header info from fileName in dataDir (optional) 
%   into a data structure V of same description as in spm12
%   NOTE: - if ANALYZE format, fileName should be *.hdr and *.img will be assumed to exist 
%         - can accomodate gziped file formats as well.
%   
%   Supported filetypes: NIFTI, ANALYZE and gzipped files
%
%   gzipFlag indicated whether to keep only archived gzipped files (lower 
%   storage cost) or allowing for ungzipped files to remain (faster run times).
%
% Jason E. Hill
% STANCE_load_header.m  updated     12 JAN 2016

% extract file extention

if nargin < 3
    gzipFlag = false;
end
if nargin < 2
    dataDir = [];
end

fileExtension = fileName(end-2:end);
gzipFileName = []; %#ok<*NASGU>
if strcmpi(fileExtension,'.gz')
    gzipFileName = fileName;
    gunzipFileName = gzipFileName(1:end-3);
    if exist(gunzipFileName,'file')
        fileName = gunzipFileName;
    else
        if nargin < 2 || isempty(dataDir)
            fileName = gunzip(fileName);
        else
            fileName = gunzip(fileName, dataDir);
        end
        fileName = char(fileName);	% convert from cell to string 
    end
end

% extract the volume header info
if nargin == 2
    filePath = fullfile(dataDir,fileName);
    V = spm_vol(filePath);
    if gzipFlag
        delete(filePath);
    end
else
    V = spm_vol(fileName);
    if gzipFlag
        delete(fileName);
    end
end

end

