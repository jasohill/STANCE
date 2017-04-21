function OUT_ts = STANCE_apply_response_function(dt, ts, response_function, param, sliceOrder)
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
%     ts:                timeseries
%     response_function: name of requested response function 
%     param:             structure of response function parameters
%     

%%
if nargin < 5
    sliceOrder = 'SD';
end
if nargin < 4
    param = [];
end
if nargin < 3
    response_function = 'canonical';
end

[Nt,ne] = size(ts.Data);
T_scan = Nt*dt;

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
    case 'ballon'
        display('o Applying the balloon response to the time-series.');         
    case {'RRF','respiratory'}
        display('o Applying the respiratory response function to time-series.'); 
    case {'CRF','cardiac'}
        display('o Applying the cardiac response function to time-series.');         
end
for i = 1:ne
    switch response_function
        case {'gamma','single','Gamma','canonical', 'double','triple','logit'}
            ts_conv_Data = STANCE_convolve_response_function(dt,ts.Data(:,ne),response_function,param,sliceOrder,30);
        case {'CRF','cardiac','RRF','respiratory'}
            ts_conv_Data = STANCE_convolve_response_function(dt,ts.Data(:,ne),response_function,param,sliceOrder,50);
        case {'Balloon','balloon'}
            ts_conv_Data = balloon(ts.Data(:,ne),dt, T_scan);
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
