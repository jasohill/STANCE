function ts_conv_Data = STANCE_convolve_response_function(dt, ts_Data, response_function, param, sliceOrder, RF_time)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%
% Generates a time series ts_conv_Data by convolving the requested response_function
% with ts. 
% The cardiac and respiratory response functions can also be called with this tool.
%
%   author:         Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
%   updated:        31 MAR 2017
%
% Arguments
%     response_function: name of requested response function 
%     param:             structure of response function parameters
%     

if nargin < 5
    sliceOrder = 'SD';
end
if nargin < 4
    param = [];
end
if nargin < 3
    response_function = 'canonical';
end
if nargin < 6
    switch response_function
            case {'gamma','single','Gamma','canonical', 'double','triple','logit'}
                RF_time = 30;
            case {'CRF','RRF','cardiac','respiratory'}
                RF_time = 50;            
    end
end
    switch response_function
            case {'gamma','single','Gamma'}
                response_function = 'singleGammaHRF';
            case {'canonical','double'}
                response_function = 'doubleGammaHRF';
            case 'triple'
                response_function = 'tripleGammaHRF';
            case 'logit'
                response_function = 'tripleLogitHRF';                
            case {'CRF','cardiac'}
                response_function = 'cardiacRF';
            case {'RRF','respiratory'}
                response_function = 'respiratoryRF';
    end
    
    Nt = length(ts_Data);

    fh = str2func(response_function);

    RF = fh((1:RF_time),'verbose',false,'param',param);
    if isvector(RF)
        RF = fh((dt:dt:RF_time),'verbose',false,'param',param);
        ts_conv_Data = conv(RF,ts_Data);  
        % impose time length
        ts_conv_Data = ts_conv_Data(1:Nt); 
        % normalize response
        ts_conv_Data = ts_conv_Data./max(ts_conv_Data);
    else % assume 4D HRF
        Nslices = size(RF,3);
        sliceTiming = make_slice_timing(sliceOrder,Nslices);
        N_ADCs = length(sliceTiming); 
        TR = N_ADCs*dt;
        Nvols = floor(Nt/N_ADCs);
        RF = fh((TR:TR:RF_time),'verbose',false,'param',param);
        ts_conv_Data = zeros(size(RF,1),size(RF,2),size(RF,3),Nvols);
        hw = waitbar(0,['Applying ', response_function,' ...']);
        for z = 1:N_ADCs
            sliceIndeces = sliceTiming(:,z);
            multiplicity = length(sliceIndeces);
            timeIndeces = z:N_ADCs:Nt;
            for s = 1:multiplicity
                sliceIndex = sliceIndeces(s); 
                RFslice = squeeze(RF(:,:,sliceIndex,:));
                ts_conv_slice = squeeze(shiftdim(convn(shiftdim(RFslice,2),ts_Data(timeIndeces)),1));
                % impose time length
                ts_conv_slice = ts_conv_slice(:,:,1:Nvols);
                % normalize response
                ts_conv_slice_max = squeeze(max(shiftdim(ts_conv_slice,2)));
                ts_conv_slice_max( ts_conv_slice_max == 0) = 1;
                for t = 1:Nvols
                    ts_conv_slice(:,:,t) = ts_conv_slice(:,:,t)./ts_conv_slice_max;
                end
                ts_conv_Data(:,:,sliceIndex,:) = ts_conv_slice;
            end
% Note: the following code solution is too demanding on storage space
%            ts_conv_Data = squeeze(shiftdim(convn(shiftdim(RF,2),ts_Data),2));
%            ts_conv_Data = ts_conv_Data(:,:,:,1:Nt);
            waitbar(z/N_ADCs,hw)
        end
    close(hw)
    end 
end