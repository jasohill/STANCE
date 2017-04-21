function ts = STANCE_event_design(dt, event_onset, event_duration, T_scan)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%
% Generates a times series of an event design
% for the specified sample time and event time onsets and durations, given
% a starting preference (default has time series start in the off position).
% 
%   Author:  Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
%   Updated: 8 MAR 2017
%
% Arguments
%     dt:    the sample time of the design [s].
% NOTE: the following if type double are in units of [s]
%                 ... if type intX/uintX are in units of dt.
% -----------------------------------------------------------------------
%     event_onset: the time onset of an event due to a single condition [s].
%     event_duration: the time interval of an event due to a single condition [s].
%     T_scan: the total time of the scan [s]. (optional)
% -----------------------------------------------------------------------
%%

% variable argument handling

if nargin < 5
    startFlag = 1;
end
if isempty(startFlag)
    startFlag = 1;
end 
if nargin < 4
    T_scan = event_onset(end) + event_duration(end);
end

if isnumeric(event_onset)
    % do nothing
elseif iscell(event_onset)
    event_onset = cell2mat(event_onset);
else
    error('The event onset times were given in an unsupported format (not cell or vector)');
end

if isnumeric(event_duration)
    % do nothing
elseif iscell(event_duration)
    event_onset = cell2mat(event_duration);
else
    error('The event duration was given in an unsupported format (not cell or vector)');
end

if isa(event_onset(1),'integer')
     event_onset = double(event_onset*dt);
end
if isa(event_duration(1),'integer')
     event_duration = double(event_duration*dt);
end
if isa(T_scan,'integer')
     T_scan = double(T_scan*dt);
end
Nevents = length(event_onset);

if (length(event_duration) == 1) && Nevents>1
    event_duration = event_duration*ones(1,Nevents);
end

%% design timeseries

% quality check
for n = 1:(Nevents-1)
    if event_onset(n+1)-event_onset(n) < event_duration(n)
        event_duration(n) = event_onset(n+1)-event_onset(n);
        warning('Event duration longer than event timing interval! ... Correcting.')
    end
end
if (T_scan - (event_onset(Nevents) + event_duration(Nevents))) < 0
    event_duration(Nevents) = T_scan-event_onset(Nevents);
    warning('Event duration longer than event timing interval! ... Correcting.')       
end

if event_onset(1) == 0.0
   startFlag = 0;  
   event_OFF_blocks = zeros(1,Nevents);
   event_OFF_blocks(1:Nevents-1) = event_onset(2:Nevents)-event_onset(1:Nevents-1) - event_duration(1:Nevents-1);
   event_OFF_blocks(Nevents)     = T_scan - (event_onset(Nevents) + event_duration(Nevents)); 
else
   startFlag = 1; 
   event_OFF_blocks = zeros(1,Nevents+1);   
   event_OFF_blocks(1)         = event_onset(1);
   event_OFF_blocks(2:Nevents) = event_onset(2:Nevents)-event_onset(1:Nevents-1) - event_duration(1:Nevents-1);
   event_OFF_blocks(Nevents+1)    = T_scan - (event_onset(Nevents) + event_duration(Nevents));   
end

ts = STANCE_design_timeseries(dt, event_duration, event_OFF_blocks, startFlag);


