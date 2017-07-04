function OUT_ts = STANCE_apply_response_function(dt, ts, response_function, param, sliceOrder, model, ratio)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%
% Generates a time series Out_ts by applying the requested response_function to ts.
% NOTE: includes the balloon model response; all others are convolutional
% The cardiac and respiratory response functions can also be called with this tool.
%
%   author:         Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
%   updated:        25 MAR 2016
%
% Arguments
%     dt:                the sample time of the timeseries [s].
%     ts:                timeseries, e.g. the experimental design boxcar function
%     response_function: name of requested response function 
%     param:             structure of response function parameters
%     sliceOrder:        the slice order string code (default = 'SD'
%                        denoting sequential descending order)
%     model:             Choices: 'sustained' (default) or 'on-off' which
%                        is a transient response depending on the derivative of ts
%     ratio:             The ratio of the onset to the offset response amplitudes
%

%%
if nargin < 7
    ratio = 1;
end
if nargin < 6
    model = 'sustained';
end
if isempty(model)
    model = 'sustained';
end
if nargin < 5
    sliceOrder = 'SD';
end
if isempty(sliceOrder)
    sliceOrder = 'SD';
end
if nargin < 4
    param = [];
end
if nargin < 3
    response_function = 'canonical';
end
if isempty(response_function)
    response_function = 'canonical';
end
    
[Nt,ne] = size(ts.Data);
T_scan = Nt*dt;

if strcmp(model,'on-off')
    % the transient model of when the experimental design turns on and off
    ts_data = ts.Data;
    for i = 1:ne
        ts_data(2:Nt,i) =  ts_data(2:Nt,i) - ts_data(1:Nt-1,i); 
    end
    if ratio < 1
        ts_data(ts_data(:,i)>0,i) = ratio*ts_data(ts_data(:,i)>0,i);  
    elseif ratio > 1 
        ts_data(ts_data(:,i)<0,i) = (1.0/ratio)*ts_data(ts_data(:,i)<0,i);         
    end
    ts.Data = abs(ts_data);
end    
    
%% generate the response output to the input time-series
switch response_function
    case {'gamma','single','Gamma'}
        display('o Applying the Gamma distribution response function to time-series.');  
    case 'canonical'
        display('o Applying the canonical hemodynamic response function to time-series.');
    case 'double'
        display('o Applying the double Gamma response function to time-series.');  
    case 'triple'
        display('o Applying the triple Gamma response function to time-series.');    
    case 'logit'
        display('o Applying the triple logit response function to time-series.'); 
    case 'balloon'
        display('o Applying the balloon response to the time-series.');         
    case {'RRF','respiratory'}
        display('o Applying the respiratory response function to time-series.'); 
    case {'CRF','cardiac'}
        display('o Applying the cardiac response function to time-series.');         
end
for i = 1:ne
    switch response_function
        case {'gamma','single','Gamma','canonical', 'double','triple','logit'}
            ts_conv_Data = STANCE_convolve_response_function(dt,ts.Data(:,i),response_function,param,sliceOrder,30);
        case {'CRF','cardiac','RRF','respiratory'}
            ts_conv_Data = STANCE_convolve_response_function(dt,ts.Data(:,i),response_function,param,sliceOrder,50);
        case {'Balloon','balloon'}
            ts_conv_Data = balloon(ts.Data(:,i),dt, T_scan);
            % no support for variability yet
    end
    
    
    if isvector(squeeze(ts_conv_Data))
        ts_conv_Data = ts_conv_Data(1:Nt);
        OUT_ts_data(:,i) = ts_conv_Data;  %#ok<*AGROW>
        DT = dt;
    elseif ndims(squeeze(ts_conv_Data)) == 4
        OUT_ts_data(:,:,:,:,i) = ts_conv_Data;
        Nvols = size(ts_conv_Data,4);
        Nslices = size(ts_conv_Data,3);
        sliceTiming = make_slice_timing(sliceOrder,Nslices);
        N_ADCs = length(sliceTiming); 
        TR = N_ADCs*dt;
        TR = N_ADCs*dt;        
        T_scan = Nvols*TR;
        DT = TR;
    end
    
end
OUT_ts_time = (DT:DT:T_scan)';
OUT_ts = timeseries(OUT_ts_data,OUT_ts_time);
OUT_ts.TimeInfo.Units = 'seconds';
OUT_ts = setuniformtime(OUT_ts,'StartTime',DT);
OUT_ts = setuniformtime(OUT_ts,'Interval',DT);
