function ts = make_impulse_timeseries(dt,Nt,TI,TI_sigma,seed)
% Generates an impulse time-series characterized by 
% the sampling time dt over a total of Nt time points,
% and average interval time TI with variance varTI.
% NOTE: the resulting time-series will have a Gaussian spectrum
% 
%   Author:  Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
%   Updated: 8 FEB 2017

% handle variable input arguments
if nargin<5
    % do nothing, carry on
else
    if ~isempty(seed)
        rng(seed,'twister');
    end
end

TT = dt*Nt;
Nk = ceil(TT*TI);
Nk = ceil(Nk*(1.1 + 2*TI_sigma)); % initial estimation of the number of intervals
TIs = TI_sigma.*randn(Nk,1) + TI;
TIts = cumsum(TIs);
lastIndex = find(TIts>TT,1);
Nk = lastIndex-1;
TIts = TIts(1:Nk);
TItis = round(TIts./dt);

%% generate impulse time-series
ts_data = zeros(Nt,1);
for i = 1:Nk
    ts_data(TItis(i)) = 1;
end
ts_time = dt:dt:(Nt*dt);
ts = timeseries(ts_data,ts_time);
ts.TimeInfo.Units = 'seconds';
ts = setuniformtime(ts,'StartTime',dt);
ts = setuniformtime(ts,'Interval',dt);

end