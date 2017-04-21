function rho_wins = find_sliding_correlation(ts,CMref,win,shift)
% Computes the Pearson's correlation coefficients for windows of size win 
% moved across ts by shift; as compared with the reference correlation 
% matrix CMref.
%
% Author: Jason E. Hill, Ph.D.
% Date: 26 Jan 2017

nT = size(ts,1);
nR = size(ts,2);

nW = 1 + floor((nT-win)/shift);

rho_wins = zeros(nW,1);

for i = 1:nW

windowed_ts = ts(shift*(i-1)+1:shift*(i-1)+win,:);

CMwin = corr(windowed_ts);

[CMrho,~] = Pearsons_correlation_coefficient(CMref,CMwin);

rho_wins(i) = CMrho;

end

end

