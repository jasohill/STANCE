function [activation_map,Y_GM] = STANCE_GM_mask(activation_map,actMNIgmVolume,This_sss,morphLimit)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% Returns the activation map masked with the gray matter tissue priors of the subject 
%
% Jason E. Hill
% STANCE_find_GM_volume.m      updated     26 FEB 2016

if nargin < 4
    morphLimit = 10;
end

if ~exist('STANCE.mat','file')
    if ~exist('../STANCE.mat','file')
       load('../../STANCE.mat'); 
    else
       load('../STANCE.mat');
    end
else
    load('STANCE.mat');
end

if nargin < 3
    This_sss = Now_sss;
end


filename_subject_GM = [STANCE_genpath(This_sss,2),'/gray_matter_mask.nii'];

[~,Y_GM] = STANCE_load_volume(filename_subject_GM);

gmVolume = sum(Y_GM(:));


activation_map0 = Y_GM.*activation_map;
% renormalize masked activation map
maxPrior = max(activation_map0(:));
activation_map0 = (activation_map0/maxPrior)...
                   .*(70./(maxPrior*70+(1-maxPrior)*60)); % partial volume adjustment: GM T2* ~70; WM/CSF T2* ~60 
ActGMvolume0 = activation_map0/gmVolume;

activation_map = activation_map0;
if sum(activation_map0(:)) > 0 
    stopFlag = false;
else
    stopFlag = true;
end
for n = 1:morphLimit
   w = 2*n + 1;
   if ~stopFlag

      activation_map_erode = imerode(activation_map0,zeros(w,w,w));     
      ActGMvolumeN = sum(activation_map_erode(:))/gmVolume;
      if (ActGMvolumeN < actMNIgmVolume) && (actMNIgmVolume < ActGMvolume0) 
          stopFlag = true;
%          -n
          if (ActGMvolume0 - actMNIgmVolume) > (ActGMvolumeN - actMNIgmVolume)    
              activation_map = activation_map_erode;
          else
              activation_map = activation_map0;              
          end
      end
      activation_map_dilate = imdilate(activation_map0,zeros(w,w,w));      
      ActGMvolumeN = sum(activation_map_dilate(:))/gmVolume;
      if (ActGMvolumeN > actMNIgmVolume ) && (actMNIgmVolume > ActGMvolume0) 
          stopFlag = true;
%          n
          if (actMNIgmVolume - ActGMvolume0) > (actMNIgmVolume - ActGMvolumeN)  
              activation_map = activation_map_dilate;
          else
              activation_map = activation_map0;              
          end
      end      
   end
end
% renormalize activation map
maxPrior = max(activation_map(:));
activation_map0 = (activation_map/maxPrior)...
                   .*(70./(maxPrior*70+(1-maxPrior)*60)); % partial volume adjustment: GM T2* ~70; WM/CSF T2* ~60 
end