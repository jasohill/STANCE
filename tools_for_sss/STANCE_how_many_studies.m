function [Nbr_studies,list] = STANCE_how_many_studies
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% Returns the number of studies
%
% Jason E. Hill
% STANCE_how_many_studies.m      updated     2 APR 2017

if ~exist('STANCE.mat','file')
    load('../STANCE.mat');
else
    load('STANCE.mat');
end

    filepath = [STANCEroot,'/fMRI'];

    if isdir(filepath)
        oldDir = cd(filepath);    
        list = ls('study*');
        Nbr_studies = size(list,1);
        cd(oldDir);
    else
        Nbr_studies = 0;
        list = [];
    end
end