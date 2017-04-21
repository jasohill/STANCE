function sliceTiming = make_slice_timing(sliceOrder,N_slices,acceleration)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%
% Generates the slice timing index vector for an MRI scan.
%
% Author: Dr. Jason E. Hill (TTU)
%
% INPUTS:
% ---------
% sliceOrder: (string specifying the slice order
%              'S_'  - Sequential ordering
%              'I_'  - Interleaved ordering
%              '_A'  - Ascending ordering
%              '_D'  - Descending ordering
%              '__2' - start with even numbers
% Nslices: (integer) 

if nargin < 3
    acceleration = 1;
end

    switch sliceOrder
        case 'SA'
            sliceTiming = (1:N_slices);        
        case 'SD'
            sliceTiming = (N_slices:-1:1);
        case 'IA'
            sliceTiming = [(1:2:N_slices) (2:2:N_slices)];            
        case 'ID'
            if mod(N_slices,2)
                % odd number of slices
                 sliceTiming = [(N_slices:-2:1) ((N_slices-1):-2:2)]; 
            else
                 sliceTiming = [((N_slices-1):-2:1) (N_slices:-2:2)];                 
            end
        case 'IA2'
            sliceTiming = [(2:2:N_slices) (1:2:N_slices)];             
        case 'ID2'
            if mod(Nslice,2)
                % odd number of slices
                 sliceTiming = [((N_slices-1):-2:2) (N_slices:-2:1)]; 
            else
                 sliceTiming = [((N_slices-1):-2:2) (N_slices:-2:1)];                 
            end            
    end

    if acceleration > 1
        sliceTiming1 = sliceTiming;
        N_ADCs = ceil(N_slices/acceleration);
        N_cols = ceil(N_slices/N_ADCs);
        if N_ADCs*N_cols > N_slices
           sliceTiming1(N_slices+1:N_ADCs*N_cols) = NaN;
        end
        sliceTiming2 = zeros(N_ADCs,N_cols);
        for c = 1:N_cols
            sliceTiming2(:,c) = sliceTiming1((1:N_ADCs)+(c-1)*N_ADCs);
        end
        sliceTiming = sliceTiming2';
    end
    
end

