function ts = make_drift_timeseries(dt,Nt,period,random_phase,order,seed)
% Generates a drift timeseries
%   author:         Dr. Jason E. Hill (post-doc fellow with CNT)
%   updated:        21 NOV 2016
%   inspiration:    somewhat by neuRosim's lowfreq() R function [Marijke Welvaert <marijke.welvaert@gmail.com>]
%   Arguments
%   dt:             Sample time [s].
%   Nt:             Number of time points.
%   period:         Characteristic period of drift oscillations [s].
%   random_phase:   Boolean flag to indicate the generation of random phases.
%   order:          The order of random polynomials to add to drift (0 =
%   none)
%   seed:           Seed for the random number generator (for reproducibility)
%   ts  :           Output is a zero mean, unit variance drift curve

% handle variable inputs
if nargin < 6
    seed = [];
end
if ~isempty(seed)
    rng(seed,'twister')
end
if nargin < 5
    order = 0;
end
if isempty(order)
    order = 0;
end 
if nargin < 4
    random_phase = false;
end
if isempty( random_phase)
    random_phase = false;
end 
if nargin < 3
    period = 128;
end
if isempty(period)
    period = inf;
end 


TT = dt*Nt;

ts_data = zeros(1,Nt);

ti = 0:(Nt-1);
if period == Inf
    % do nothing
else
    % compute number of basis functions
    Nbasis = floor((2*TT)/period + 1);
    if Nbasis < 3
        Nbasis = 3;
        warning('May mis-model drift! A longer scan time OR a lower period may work better.');
    end
    for k = 2:Nbasis
        if random_phase
            phi = 2*pi*rand;
        else
            phi = 0;
        end
        ts_data = ts_data + cos(pi*(2*ti+1)*(k-1)/(2*Nt)+phi);
    end
end
ts_span = max(ts_data) - min(ts_data);

order = round(abs(order));  
if order == 0
    %do nothing
else
    for p = 1:order
        a_p = max([ts_span,1.0])*(2.0*rand-1.0);   % coefficient in (-1,+1)
        lb  = -1.0*rand;        % lower bound in (-1,0)
        ub  =  1.0*rand;        % upper bound in (0,+1)
        ts_data = ts_data + a_p*((ub-lb)*ti/(Nt-1)+lb).^p;
    end
end

ts_data = (ts_data-mean(ts_data))/std(ts_data);

ts_time = dt:dt:(Nt*dt);
ts = timeseries(ts_data,ts_time);
ts.TimeInfo.Units = 'seconds';
ts = setuniformtime(ts,'StartTime',dt);
ts = setuniformtime(ts,'Interval',dt);

end

