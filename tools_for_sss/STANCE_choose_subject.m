function [V,Y] = STANCE_choose_subject(s,mode,gzipflag)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% Select the subject brain to use for STANCE.
%
% Jason E. Hill
% STANCE_choose_subject.m      updated     21 JUN 2016

% handle variable arguments
if nargin < 3
    gzipflag = false;  % this allows for the creation of an unzipped file
end

if nargin < 2
    mode = 'fuzzy';
end

% handle file structure
addpath(genpath(pwd))

if ~exist('STANCE.mat','file')
    load('..\STANCE.mat');
else
    load('STANCE.mat');
end

% handle mode
if strcmpi(mode,'fuzzy')
    postfix = subjectsFuzzyPostfix;
elseif strcmpi(mode,'T1')
    postfix = subjectsT1postfix;
else
    warning('Unrecognized mode specified')
    postfix = [mode,'.nii'];
end

% build filepath
if s < 4
    subject  = [subjectsBase,'0', num2str(subject_labels(s)),'/subject0', num2str(subject_labels(s))];  
else
    subject  = [subjectsBase, num2str(subject_labels(s)),'/subject', num2str(subject_labels(s))];          
end           
fileNameNII  = [subject,postfix];   

% load subject brain

[V,Y] = STANCE_load_volume(fileNameNII,[],gzipflag);

end

