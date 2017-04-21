function combo_map = STANCE_combine_maps(operation, A, varargin)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% generate activation region for specified dimensions
% 
%  fuzzy logical operation followed by a list of maps: A, B, etc.
%  NOTE: if only A then taken as an array of maps already constructed,
%  unless operation = 'NOT'
%
%  operation:  'NOT'        (1.0-A)
%              'OR'         max(A,B,...)
%              'AND'        min(A,B,...)
%              'XOR'        max(A.*(1-B)...,B.*(1-A).., ...)
%              'NAND'/'MINUS'/'SUBTRACT'
%                          A.*(1-B).*...
%
% Jason E. Hill
% STANCE_combine_maps.m      updated     11 SEPT 2016

sizeA = size(squeeze(A));
if nargin < 3
   maps = A; 
else
    nVarargs = length(varargin);
    maps(1,:) = A(:);
    for m = 1:nVarargs
        B = varargin{m};
        maps(1+m,:) = B(:);
    end
end
Nmaps = size(maps,1);
    
switch upper(operation)
    % only one input:
    case 'NOT'
        combo_map = (1.0-A);
    % any number of inputs:
    case 'OR'
        OR_map = max(maps);
        combo_map = reshape(OR_map,sizeA);
    case 'AND'    
        AND_map = min(maps);    
        combo_map = reshape(AND_map,sizeA);        
    case 'XOR'
        for m = 1:Nmaps
            for n = 1:Nmaps
                if n ~= m 
                    maps(m,:) = maps(m,:).*(1-maps(n,:));
                end
            end
        end
        XOR_map = max(maps);
        combo_map = reshape(XOR_map,sizeA);  
    case {'NAND','MINUS','SUBTRACT'}
        NAND_map = maps(1,:);
        for n = 2:Nmaps
            NAND_map = NAND_map.*(1-maps(n,:));
        end        
        combo_map = reshape(NAND_map,sizeA);            
end

% renormalize the combination map
maxMap = max(combo_map(:));
if maxMap > 0;
    combo_map = combo_map/maxMap;
end


end
