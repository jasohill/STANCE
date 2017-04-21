function map = STANCE_parse_combine(task)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% combine the activation maps specified in task according to
% the operations and options in task.combine
%   
% Jason E. Hill
% STANCE_parse_combine.m      updated     22 SEP 2016

dimensions = size(task.activation(1).map);

if ~isfield(task,'combine')
    if length(task.activation) == 1
        map = task.activation.map;
    else
        % DEFAULT: combine all components with fuzzy logic OR
        map = STANCE_combine_maps('OR',task.activation(:).map);
    end
else
    for c = 1:length(task.combine)
         operation = task.combine{c}{1};
         option = task.combine{c}{2};
         if strcmp(operation,'mask')
             if strcmp(option,'L')
                 % restrict activation to left side
                 map(ceil(0.5*dimensions(1)):end,:,:) = 0;
             elseif strcmp(option,'R')
                 % restrict activation to right side
                 map(1:floor(0.5*dimensions(1)),:,:) = 0;                 
             else
                 warning('Unknown mask option.');
             end
         else
             % assume fuzzy logic operation
             if strcmp(option,'all')           
                 map = STANCE_combine_maps(operation,task.activation(:).map);
             elseif strcmp(option,'flip') 
                 map = STANCE_combine_maps(operation,map,flipud(map));             
             else
                 if length(option) == 1
                     map = STANCE_combine_maps(operation,task.activation(option(1)).map);                   
                 else
                     map = task.activation(option(1).map);
                     for a = 2:length(option)
                         map = STANCE_combine_maps(operation,map,task.activation(option(a)).map);                   
                     end
                 end
             end
         end
    end
end

end