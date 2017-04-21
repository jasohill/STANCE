function s = STANCE_random_subject(sex)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% Select a random subject brain to use for STANCE.
%
% Jason E. Hill
% STANCE_new_study.m      updated     25 FEB 2016

if nargin < 1
    sex = 'any';  % this allows for the random selection of any subject
end

if ~exist('STANCE.mat','file')
    load('..\STANCE.mat');
else
    load('STANCE.mat');
end

% select subject index s
if strcmpi(sex,'male')
    s = male_labels(randi(length(male_labels)));
elseif strcmpi(sex,'female')
    s = female_labels(randi(length(female_labels))); 
else
    s = randi(Nbr_subject_labels);
end    

end

