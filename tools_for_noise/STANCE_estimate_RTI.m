function [RTI_ave,RTI_sigma] = STANCE_estimate_RTI(y,TR)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE)
%
% Estimates the average respiratory time interval from the spectral features
% of a motion time series (usually y-direction). 
% TR is the repetition time [s]
%
%   author:         Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
%   updated:        5 APR 2017

Y = fft(y,1001);
Pyy = Y.*conj(Y)/1001;
fs = (0:1000)/(TR*1001);
max_fs = fs(500);
if max_fs > 0.275;
    hi_idx = find(fs>0.275,1);
else
    hi_idx = 500;
end
mid_idx = find(fs>0.225,1);
lo_idx = find(fs>(2*0.225-min([0.275 max_fs])),1);
fs_Resp_win = fs(lo_idx:hi_idx);
Pyy_Resp_win = Pyy(lo_idx:hi_idx);
max_Resp_freq = max(Pyy_Resp_win);
max_Resp_idx = find(Pyy_Resp_win == max_Resp_freq,1);
RTIs = 1./fs_Resp_win(Pyy_Resp_win>0.75*max_Resp_freq)';
% compute average time interval from weighted average
RTI_ave = sum(RTIs.*Pyy_Resp_win(Pyy_Resp_win>0.75*max_Resp_freq))/sum(Pyy_Resp_win(Pyy_Resp_win>0.75*max_Resp_freq));
if length(RTIs)>2
    RTI_sigma = std(RTIs);
else
    RTI_sigma = std(1./fs_Resp_win(Pyy(lo_idx:hi_idx)>0.5*max_Resp_freq));
end

end

