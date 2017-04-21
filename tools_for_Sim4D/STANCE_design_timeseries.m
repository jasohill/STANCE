function ts = STANCE_design_timeseries(dt, event_ON_blocks, event_OFF_blocks, startFlag)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%
% Generates a time series ts with boxcar function transitions
% for the specified sample time and event on and off time intervals, given
% a starting preference (default has time series start in the off position).
% 
%   Author:        Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
%
% Arguments
%     dt:              the sample time of the design [s].
%
% NOTE: the following if type double are in units of [s]
%                 ... if type intX/uintX are in units of dt.
% -----------------------------------------------------------------------
%     event_blocks:   the time intervals for events.
%     off_blocks:     the time intervals between events.
% -----------------------------------------------------------------------
%     startFlag:       whether to start with the off timeseries (default = true).
%%

% variable argument handling

if nargin < 4
    startFlag = 1;
end
if isempty(startFlag)
    startFlag = 1;
end 

% error checking
if isnumeric(event_ON_blocks)
    Non = length(event_ON_blocks);
elseif iscell(event_ON_blocks)
    event_ON_blocks = cell2mat(event_ON_blocks);
    Non = length(event_ON_blocks);
else
    error('Event ON timeseries input was given in an unsupported format (not cell or vector)');
end

if isnumeric(event_OFF_blocks)
    Noff = length(event_OFF_blocks);
elseif iscell(event_OFF_blocks)
    event_OFF_blocks = cell2mat(event_OFF_blocks);
    Noff = length(event_OFF_blocks);
else
    error('Event OFF timeseries input was given in an unsupported format (not cell or vector)');
end

if  Noff == 1
    if startFlag
        event_OFF_blocks = event_OFF_blocks.*ones(1,length(Non+1));
    else
        event_OFF_blocks = event_OFF_blocks.*ones(1,length(Non));        
    end
elseif ~(Noff == Non || Noff == (Non+1)) && startFlag
    error('There is a mismatch between the number of ON and OFF blocks.');
elseif ~(Noff == (Non-1)  || Noff == Non) && ~startFlag
    error('There is a mismatch between the number of ON and OFF blocks.');
end

if isa(event_ON_blocks,'integer')
     event_ON_blocks = double(event_ON_blocks*dt);
end
if isa(event_OFF_blocks,'integer')
     event_OFF_blocks = double(event_OFF_blocks*dt);
end

%% initializations

Nt = round((sum(event_ON_blocks) + sum(event_OFF_blocks))/dt);
ts_data = zeros(Nt,1);
event_onsets = ones(Non,1);

if startFlag
    event_onsets(1) = round(event_OFF_blocks(1)/dt);
    if event_onsets(1) < 1
        event_onsets(1) = 1;
    end
    if event_onsets(end) > Nt
        event_onsets(end) = Nt;
    end     
    for i = 2:Non
        event_onsets(i) = round((sum(event_ON_blocks(1:i-1)) + sum(event_OFF_blocks(1:i)))/dt);
        if event_onsets(i) < 1
            event_onsets(i) = 1;
        end
        if event_onsets(i) > Nt
            event_onsets(i) = Nt;
        end        
    end
else
    event_onsets(1) = 1;
    for i = 2:Noff
        event_onsets(i) = round((sum(event_ON_blocks(1:i-1)) + sum(event_OFF_blocks(1:i-1)))/dt);
        if event_onsets(i) < 1
            event_onsets(i) = 1;
        end  
        if event_onsets(i) > Nt
            event_onsets(i) = Nt;
        end         
    end    
end

%% design timeseries

for i = 1:length(event_ON_blocks)
    ts_data(event_onsets(i):(event_onsets(i) + round(event_ON_blocks(i)/dt))) = 1;
end

% correct incorrect length due to any rounding
ts_data = ts_data(1:Nt);
ts_time = dt:dt:(Nt*dt);
ts = timeseries(ts_data,ts_time);
ts.TimeInfo.Units = 'seconds';
ts = setuniformtime(ts,'StartTime',dt);
ts = setuniformtime(ts,'Interval',dt);
