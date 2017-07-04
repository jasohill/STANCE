function ts = make_ARMA_timeseries(dt,Nt,sigma,p,q,constant)
% Generates the ARMA model time correlations to a timeseries
%   author:         Dr. Jason E. Hill (post-doc fellow with CNT)
%   updated:        21 NOV 2016
%   inspiration:    somewhat by neuRosim's temporalnoise() R function [Marijke Welvaert <marijke.welvaert@gmail.com>]
%   Arguments
%   dt:        Sample time.
%   Nt:        Number of time points.
%   sigma:     The standard deviation of the signal.
%   for variable options:   'AR'         AR order followed by value (default = 0.2)
%                           'MA'         MA order followed by value (default = 0.0)
%                           'Constant'   ARMA constant followed by value (default = 0.0)
%


%% handle variable argument options

if nargin < 6
    constant = 0.0;
end
if isempty(constant)
    constant = 0.0;
end 

if nargin < 5
    q = 0.0;
end
if isempty(q)
    q = 0.0;
end 

if nargin < 4
    p = 0.2;
end
if isempty(p)
    p = 0.2;
end 

%% Construct ARMA model and timeseries

ARMA_model = arima('AR',p,'MA',q,'Variance',sigma.^2,'Constant',0);

[ts_data,~] = simulate(ARMA_model,Nt);
 
ts_time = dt:dt:(Nt*dt);
ts = timeseries(ts_data,ts_time);
ts.TimeInfo.Units = 'seconds';
ts = setuniformtime(ts,'StartTime',dt);
ts = setuniformtime(ts,'Interval',dt);

end

