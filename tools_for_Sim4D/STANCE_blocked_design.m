function ts = STANCE_blocked_design(dt, start_lag, block, ISI, T_scan, startFlag)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%
% Generates a times series of a block design
% for the specified sample time and on and off time durations, given
% a starting preference (default has time series start in the off position).
% 
%   Author:  Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
%   Updated: 25 OCT 2016
%
% Arguments
%     dt:    the sample time of the design [s].
% NOTE: the following if type double are in units of [s]
%                 ... if type intX/uintX are in units of dt.
% -----------------------------------------------------------------------
%     start_lag: zero padding at the beginning of the time-series to
%     accomodate other conditions. If type = double [s]; type = integer [# of dt].
%     block: the time interval of an event due to a single condition. If type = double [s]; type = integer [# of dt].
%     ISI:   the inter-stimulus interval. If type = double [s]; type = integer [# of dt]. 
%     T_scan: the total time of the scan. If type = double [s]; type = integer [# of dt].
% -----------------------------------------------------------------------
%     startFlag:   whether to start with the off timeseries (default = true).
%%

% variable argument handling

if nargin < 6
    startFlag = 1;
end
if isempty(startFlag)
    startFlag = 1;
end 

% error checking
if isnumeric(start_lag)
    % do nothing
elseif iscell(start_lag)
    start_lag = cell2mat(start_lag);
else
    error('The start lag time was given in an unsupported format (not cell or vector)');
end
Nevents = length(start_lag);

if isnumeric(block)
    % do nothing
elseif iscell(block)
    block = cell2mat(block);
else
    error('The block (event interval) was given in an unsupported format (not cell or vector)');
end

if isnumeric(ISI)
    % do nothing
elseif iscell(ISI)
    ISI = cell2mat(ISI);
else
    error('The ISI was given in an unsupported format (not cell or vector)');
end

if isa(start_lag,'integer')
     start_lag = double(start_lag*dt);
end
if isa(block,'integer')
     block = double(blocks*dt);
end
if isa(ISI,'integer')
     ISI = double(ISI*dt);
end
if isa(T_scan,'integer')
     T_scan = double(T_scan*dt);
end

if (length(block) == 1) && Nevents>1
    block = block*ones(1,Nevents);
end
Nevents = length(block);
if (length(ISI) == 1) && Nevents>1
    ISI = ISI*ones(1,Nevents);
end

%% design timeseries

Ncycles = floor((T_scan - start_lag)./(block + ISI));

for n = 1:Nevents

    if (T_scan(n)-start_lag(n)-Ncycles(n)*(block(n) + ISI(n))) > block(n) 
        Non = Ncycles(n) + 1;
    else
        Non = Ncycles(n);
    end
    event_ON_blocks = block(n).*ones(Non,1);
    if startFlag
        event_OFF_blocks = ISI(n).*ones(Non+1,1);
        event_OFF_blocks(1) = start_lag;
    else
        event_OFF_blocks = ISI(n).*ones(Non,1);
    end
    T_pad = (T_scan(n) - start_lag(n) - Ncycles(n)*ISI(n) - Non*block(n));
    if T_pad >0
        event_OFF_blocks(end) = (T_scan(n) - start_lag(n) - Ncycles(n)*ISI(n) - Non*block(n));    
    end
    if event_OFF_blocks(end) == 0
        event_OFF_blocks = event_OFF_blocks(1:end-1);
    end
    
    if n == 1
        ts = STANCE_design_timeseries(dt, event_ON_blocks, event_OFF_blocks, startFlag);
        ts_data(:,1) = ts.Data;
    else
        tsn = STANCE_design_timeseries(dt, event_ON_blocks, event_OFF_blocks, startFlag);
        ts_data(:,n) = tsn.Data;         %#ok<AGROW>
    end
end

if Nevents > 1
    ts.Data = ts_data;
end


