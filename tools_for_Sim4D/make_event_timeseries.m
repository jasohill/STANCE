function [ts,Nk] = make_event_timeseries(ms,TI)
% Generates an impulse event time-series characterized by 
% the modulation function time-series ms and average interval time (TI) 
% according to the Intergral Pulse Frequency Modulation (IPFM) model,
% where the threshold is normalized to R = 1.
% NOTE: the resulting time-series will have a spectrum dictated by ms.
% 
%   Author:  Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
%   Updated: 8 FEB 2017

ms_Data = ms.Data;
dt = ms.TimeInfo.Increment;
Nt = length(ms_Data);

%% Generate event impulse timeseries
ts_data = zeros(1,Nt);

n = 1;
k = 1;
while n < Nt
    t_integral = 0.0;
    while (t_integral < 1.0) && (n < Nt)
        t_integral = t_integral + ms.Data(n)/TI; 
        n = n+1;
    end
    ts_data(n) = 1;
    TIts(k) = n;
    Nk = k;    
    k = k+1;
end

TItis = TIts.*dt;

ts_time = dt:dt:(Nt*dt);
%size(ts_data)
%size(ts_time)
ts = timeseries(squeeze(ts_data),ts_time);
ts.TimeInfo.Units = 'seconds';
ts = setuniformtime(ts,'StartTime',dt);
ts = setuniformtime(ts,'Interval',dt);

end