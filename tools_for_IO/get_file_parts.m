function [fileBase,fileExtension] = get_file_parts(fileName)
% Returns base of file name and file extension (unzipped)
%
% Jason E. Hill
% get_file_parts.m      updated     27 FEB 2016


fileExtension = fileName(end-2:end);
if strcmpi(fileExtension,'.gz')
    fileBase = fileName(1:end-7);
    fileExtension = fileName(end-6:end-3);
else
    fileBase = fileName(1:end-4);  
    fileExtension = fileName(end-3:end);    
end

end

