function activation_map = STANCE_make_activation_map(dimensions, origin, activation, method)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% generate activation region for specified dimensions
%   
%   inputs arguments:
%   dimensions:(vector) Dimensions of the image space
%   origin:(vector)     Stereotactic origin of the space. Give [] to specify coordinates as is   
%   ----------------------
%   activation structure:
%   ----------------------
%   a.region:(char)     Name of the region, atlas or source filename.
%
%   a.center:(vector)   Stereotactic coordinates of the activated region center relative to origin,
%                       if form is "cube" or "sphere" the coordinates represent the center of the region,
%                       if form is "mask"/"manual", the coordinates should be in matrix
%                       form, where the rows represent the voxels and the
%                       columns the x-y/x-y-z coordinates.
%                      
%   a.volume:(double)   The volume of the activated region [mm^3] or
%                       ROI/component label in atlas/file.
%                       If absent use ratio as sides/radii with shape to
%                       calculate volume
%                       Note: volume < 0 also means to use ratio vector as the
%                       length of sides/radii of the shape as specified by
%                       'shape.' In that case if ratio unspecified, abs(volume) is used to
%                       generate it.
%   
%   a.proportion:(vec)  default [1,1]/[1,1,1], the 2D/3D aspect ratio of sides/radii
%                       e.g. 'ellipsoid' [1,2,4] has a 1:2:4 aspect ratio, i.e. a = 4c, b = 2c.  
%                       NOTE: if volume not specified use as length of
%                       sides/radii along with shape to generated volume.
%
%   a.shape:(char)      The shape of the activated regions (or ROI name):
%   a.shape{1}:(cell)   -----------------------------
%                       2-D shapes (Area -> volume)
%                       -----------------------------
%                       'square':    Area = s^2         (also 'cube') 
%                       'rectangle': Area = a*b         (also 'cuboid')
%                       'circle':    Area = 4*pi*r^2    (also 'sphere')  
%                       'ellipse':   Area = 4*pi*a*b    (also 'ellipsoid')
%                       'supercircle':  n = abs(superpower) (also 'supersphere')
%                                    Area = 4*r^2*gamma(1+1/n)^2/gamma(1+2/n)
%                       'squircle':     n = 4
%                       'diamond':      n = 1
%                                    Area = (1/(4*sqrt(2*pi))r^2*gamma(1/4)^2
%                       'astroid':      n = 2/3
%                                    Area = (3*pi/8)r^2
%                       'superellipse': n = abs(volume) (also 'superellipsoid')
%                                    Area = 4*a*b*gamma(1+1/n)^2/gamma(1+2/n)
%                       ---------------------------
%                       3-D shapes (V -> volume)
%                       ---------------------------
%                       'cube':      V = s^3     
%                       'prism':     V = a^2*b
%                       'cuboid':    V = a*b*c          (also 'box')
%                       'sphere':    V = (4/3)*pi*r^3   (also 'ball')
%                       'spheroid':  V = (4/3)*pi*a^2*b
%                       'ellipsoid': V = (4/3)*pi*a*b*c
%                       'supersphere':  n = abs(superpower) (also 'superball')
%                                    V = (8/(3*n))*r^3*beta(2/n,1/n) 
%                       'squircle':     n = 4
%                       'diamond':      n = 1
%                       'astroid':      n = 2/3
%                                    V = (64/105)*r^3
%                       'superellipsoid' 
%                             case: length(superpower) = 1 
%                                       n = abs(superpower)
%                                    V = (8/(3*n))*a*b*c*beta(2/n,1/n)
%                             case: length(superpower) = 2 
%                                       n = abs(superpower(1)), r = abs(superpower(2))
%                                    V = (8/(3*r*n))*a*b*c*beta(1/r,1/r)*beta(2/n,1/n)
%                       'superegg' 
%                                       n = abs(superpower)
%                                    V = (8*pi/(6*n))*r^2*h*beta(2/n,1/n)
%                       ---------------------------------------
%                       Defined manually by user-defined mask
%                       ---------------------------------------
%                       'mask' (also 'manual')
%                       ---------------------------------------
%                       Defined by atlas ROI 
%                       ---------------------------------------
%                       'atlas' (also 'parcellation')
%
%   shape{2}:(d/v)     The exponent powers to use with 'super' shapes.
%                      abs(power) = the power n
%                      or it is a vector with absolute value [t e]
%                      The threshold for altas.
%
%   a.falloff:(double)  Decay rate between 0 and 1. 0 means no falloff,
%                       while 1 results in the fastest decay.
%
%   a.minimum:(double)  Minimum value of activation with template mask: between 0 and 1. 
%
%   a.rotation:(d/v)    If 2D, an angle value in degrees. 
%                       If 3D, the Tait-Bryan angles of alpha, beta, gamma for
%                       X,Y,Z axis rotation individually (pitch, roll, yaw) in degrees. 
%
%   method              Interpolation method. (Default = 'linear')
%
%   output:
%   act:           An array representing the activation image with specified regions
%
% Xiangyu Liu & Jason E. Hill
% STANCE_make_activation_map.m      updated     28 MAR 2017


%% optional argument handling

if nargin < 4
   method = 'linear';
end
Ndims = length(dimensions);

if isfield(activation,'region')
    region = activation.region;
else
    region = '';
end

if isfield(activation,'center')
    center = activation.center;
else
    center = [];
end
if isfield(activation,'volume')
    volume = activation.volume;
else
    volume = [];
end
if isfield(activation,'rotation')
    rotation = activation.rotation;
else
    rotation = [];
end
if isfield(activation,'proportion')
    proportion = activation.proportion;
else
    proportion = [];
end
if isfield(activation,'falloff')
    falloff = activation.falloff;
else
    falloff = [];
end
if isfield(activation,'minimum')
    minimum = activation.minimum;
else
    minimum = [];
end

if iscell(activation.shape)
    shape      = activation.shape{1};
    superpower = activation.shape{2}; 
else
    shape      = activation.shape;
end

if isempty(origin)
    if Ndims == 2
        origin = [0,0];        
    else
        origin = [0,0,0];
    end
end

if isempty(center)
    if Ndims == 2
        center = [0,0];        
    else
        center = [0,0,0];
    end
end

if isempty(shape)
    if Ndims == 2
        shape = 'circle';        
    else
        shape = 'sphere';
    end
end

if isempty(proportion)
    if isempty(volume)
        if Ndims == 2
            proportion = [1, 1];        
        else
            proportion = [1, 1, 1];
        end
    elseif volume>0
        if Ndims == 2
            proportion = [1, 1];        
        else
            proportion = [1, 1, 1];
        end
    else
        if Ndims == 2
            proportion = [abs(volume), abs(volume)];        
        else
            proportion = [abs(volume), abs(volume), abs(volume)];
        end
    end
else  % swap XY coordinate directions
    if length(proportion) > 1
        temp = proportion(1);
        proportion(1) = proportion(2);
        proportion(2) = temp; 
    end
end

if isempty(rotation)
    if Ndims == 2
        rotation = 0;
    else
        rotation = [0, 0, 0];
    end
end

if isempty(falloff)
    falloff = 0;
end

if isempty(minimum)
    minimum = 0;
end

if Ndims<2 || Ndims>3
    error('Dimensions should represent a 2D image or 3D volume');
end

if strcmp(shape, 'cube') || strcmp(shape, 'sphere') || strcmp(shape, 'square') || strcmp(shape, 'rectangle') || strcmp(shape, 'box') || strcmp(shape, 'circle') || strcmp(shape, 'ball') || strcmp(shape,'ellipse') || strcmp(shape, 'diamond') || strcmp(shape,'cuboid') || strcmp(shape,'prism') ||  strcmp(shape,'spheroid') ||  strcmp(shape,'ellipsoid') || strcmp(shape, 'superellipse') || strcmp(shape, 'squircle') || strcmp(shape, 'superellisoid') || strcmp(shape, 'astroid')  || strcmp(shape, 'superball') || strcmp(shape, 'superegg')     
    if length(center) ~= length(dimensions) || length(origin) ~= length(dimensions)
        error('Mismatch between dimensions of target space and coordinates');
    end
end

if strcmp(shape, 'square') || strcmp(shape, 'rectangle') || strcmp(shape, 'circle') || strcmp(shape, 'superellipse') 
    if Ndims == 3
        error('Mismatch between dimensions of image space and coordinates');
    end
end

if strcmp(shape, 'mask') || strcmp(shape, 'manual')
    if ismatrix(center)
        if size(center,2) ~= length(dimensions)
            error('Mismatch between dimensions of image space and coordinates');
        end
    else
        error('The coordinates should be a matrix');
    end
elseif ~strcmp(shape, 'atlas') || strcmp(shape, 'parcellation')
    shiftVector = round(center+origin);
    if abs(sum((center+origin)-shiftVector)) > 0
        % handle non-integer shifts
        residualShift = (center+origin)-shiftVector;
    else
        if Ndims == 2
            residualShift = [0, 0];
        else
            residualShift = [0, 0, 0];
        end
    end
end

% extract exponents for 'super' shapes
switch shape
    case {'supercircle','supersphere','superball','superellipse','superellipsoid'}
        n = abs(superpower(1));
        if length(superpower) == 2
            r = abs(superpower(2));
        else
            r = 1;
        end
    case {'squircle'}    
        n = 4;
        r = 4;
    case {'diamond'}            
        n = 1;
        r = 1;
    case {'astroid'}            
        n = 2/3;
        r = 2/3;
    case {'superegg'}            
        n = abs(superpower(1));
        r = 2;
end
if ~isempty(volume)
    if volume < 0 && ~strcmp(shape,'atlas') && ~strcmp(shape,'parcellation')
           volume = [];
    end
else
    if Ndims == 2
    switch shape
        case {'square','cube'}
            volume = proportion(1)^2;
            proportion = [1 1];
        case {'rectangle','cuboid','box'}
            if length(proportion)>1
                volume = proportion(1)*proportion(2);
            else
                volume = proportion^2;                
            end
        case {'circle','sphere','ball'}
            volume = 4*pi*proportion(1)^2;
            proportion = [1 1];
        case {'ellipse','ellipsoid','spheroid'}
            if length(proportion)>1
                volume = 4*pi*proportion(1)*proportion(2);
            else
                volume = 4*pi*proportion(1)^2;
            end
        case {'supercircle','supersphere','superball','squircle','astroid','diamond'}
            volume = (4*gamma(1+1/n)^2/gamma(1+2/n))*proportion(1)^2;
        case {'superellipse','superellipsoid'}
            volume = (4*gamma(1+1/n)^2/gamma(1+2/n))*proportion(1)*proportion(2);
        case {'mask','manual'}
            volume = size(center,1);
        case {'atlas','parcellation'}    
            %do nothing
        otherwise
            error('Unsupported shape, cannot compute activation volume.')
    end
    else
    switch shape
        case {'cube'}
            volume = proportion(1)^3;
            proportion = [1 1 1];
        case {'prism'}
            if length(proportion)>1            
                volume = proportion(1)^2*proportion(2);
                proportion = [1 1 proportion(2)/proportion(1)];
            else
                volume = proportion(1)^3;
            end
        case {'cuboid','box'}
            if length(proportion)>2
                volume = proportion(1)*proportion(2)*proportion(3);
            elseif size(proportion)>1
                volume = proportion(1)^2*proportion(2);                 
            else
                volume = proportion^3;                
            end
        case {'sphere','ball'}
            volume = (4*pi/3)*proportion(1)^3;
            proportion = [1 1 1];            
        case {'spheroid'}
            if length(proportion)>1
                volume = (4*pi/3)*proportion(1)^2*proportion(2); 
                proportion = [1 1 proportion(2)/proportion(1)];                
            else
                volume = (4*pi/3)*proportion(1)^3;
            end
        case {'ellipse','ellipsoid'}
            if length(proportion)>1
                volume = (4*pi/3)*proportion(1)*proportion(2)*proportion(3); 
            elseif length(proportion)>2
                volume = (4*pi/3)*proportion(1)^2*proportion(2); 
            else
                volume = (4*pi/3)*proportion(1)^3;
            end           
        case {'supersphere','superball','astroid','squircle','diamond'}
            volume = ((8/(3*n))*beta(2/n,1/n))*proportion(1)^3;
        case {'superegg'}
            volume = (8*pi/(6*n))*beta(2/n,1/n)*proportion(1)^2*proportion(2);
        case {'superellipsoid'}
            volume = (8/(3*r*n))*beta(1/r,1/r)*beta(2/n,1/n)*proportion(1)*proportion(2)*proportion(3);            
        case {'mask','manual'}
            volume = size(center,1);
        case {'atlas','parcellation','data'}    
            %do nothing
        otherwise
            error('Unsupported shape, cannot compute activation volume.')             
    end
    end
end

if ~isempty(rotation)
    rotation = rotation*pi/(180);
    if Ndims==2
        if length(rotation) ~= 1
            error('Mismatch between dimensions of image space and rotation option');
        end
        Theta = rotation;
        Rotation = [cos(Theta), sin(Theta), 0; -sin(Theta), cos(Theta), 0; 0,0,1];
    else
        if length(rotation) ~= 3
            error('Mismatch between dimensions of image space and rotation option');
        end
        TaitBryanAlpha = rotation(2);
        TaitBryanBeta  = rotation(1); 
        TaitBryanGamma = rotation(3);
        Rotation = [cos(TaitBryanBeta)*cos(TaitBryanGamma),cos(TaitBryanGamma)*sin(TaitBryanAlpha)*sin(TaitBryanBeta)-cos(TaitBryanAlpha)*sin(TaitBryanGamma),cos(TaitBryanAlpha)*cos(TaitBryanGamma)*sin(TaitBryanBeta)+sin(TaitBryanAlpha)*sin(TaitBryanGamma),0;...
            cos(TaitBryanBeta)*sin(TaitBryanGamma), cos(TaitBryanAlpha)*cos(TaitBryanGamma)+sin(TaitBryanAlpha)*sin(TaitBryanBeta)*sin(TaitBryanGamma),-cos(TaitBryanGamma)*sin(TaitBryanAlpha)+cos(TaitBryanAlpha)*sin(TaitBryanBeta)*sin(TaitBryanGamma),0;...
            -sin(TaitBryanBeta),cos(TaitBryanBeta)*sin(TaitBryanAlpha),cos(TaitBryanAlpha)*cos(TaitBryanBeta),0;...
            0,0,0,1];
        Rx = [1,0,0,0;0,cos(TaitBryanAlpha),-sin(TaitBryanAlpha),0;0,-sin(TaitBryanAlpha),cos(TaitBryanAlpha),0;0,0,0,1];
        Ry = [cos(TaitBryanBeta),0,-sin(TaitBryanBeta),0;0,1,0,0;sin(TaitBryanBeta),0,cos(TaitBryanBeta),0;0,0,0,1];
        Rz = [cos(TaitBryanGamma),sin(TaitBryanGamma),0,0;-sin(TaitBryanGamma),cos(TaitBryanGamma),0,0;0,0,1,0;0,0,0,1];
    end
end

%%
activation_map = zeros(dimensions);

% define encompassing dimensions: ~2x extent to allow for centers near edges
ratioMinDim = find(proportion == min(proportion),1,'first');
encompassingDims = 2*round(2.*(floor(proportion(ratioMinDim)*0.5.*dimensions./proportion)+0.5))+1;

if Ndims == 2
    switch shape
        
        case {'square','rectangle','cube','box'} % a square/rectangle in 2D
            halfsides = 0.5*sqrt(volume/prod(proportion)).*proportion;
            idx = (halfsides==0);
            halfsides(idx)=1;
            % first, specify a square region
            for d=1:2
               % region dimensions
               regionDims(d) = min(2*(ceil(min(halfsides)))+3,encompassingDims(d));
               % center coordinate of sphere template
               squareCenter(d) = (regionDims(d)+1)/2;
            end
            img = zeros(regionDims);             
            % next, specify an exact square with ratio 1:1 with smallest side
            disp('o Specifying square template')            
            squareHalfSide = min(halfsides) + 1;
            squareCenter = (regionDims+1)/2; % square template center 
            for i = 1:regionDims(1)
                for j = 1:regionDims(1)
                    if abs(i-squareCenter(1)+residualShift(1))<=squareHalfSide && abs(j-squareCenter(2)+residualShift(2))<=squareHalfSide
                        img(i,j) = 1;
                    end
                end
            end
            disp('o Performing affine transformation.')            
            ratioMinDim = find(halfsides == min(halfsides),1,'first'); 
            Scale = [proportion(1)/proportion(ratioMinDim),0,0;0,proportion(2)/proportion(ratioMinDim),0;0,0,1];
            M = Scale'*Rotation';
            tshape = affine2d(M);
            RI = imref2d(size(img),1,1);
            newImg = imwarp(uint8(255*img),RI,tshape,method);
            newDims = size(newImg);
            padDims = ceil(0.5*(dimensions-newDims));         
            padDims(padDims<0) = 0;
            newImg = padarray(newImg,padDims);
            newDims = size(newImg);
            newDims = 2.*(ceil(0.5.*newDims));            
            disp('o Building activation map.') 
            newImg = circshift(newImg,shiftVector-newDims/2);            
            center = origin + center;
            for i = 1:dimensions(1)
                for j = 1:dimensions(2)
                    if newImg(i,j) > 0
                        if falloff ~= 0
                            activation_map(i,j) = (double(newImg(i,j))/255.0)*((1-minimum)*exp(-(((i-center(1))^2+(j-center(2))^2))*falloff)+minimum);
                        else
                            activation_map(i,j) = double(newImg(i,j))/255.0;
                        end
                    end
                end
            end
            
        case {'circle','ellipse','sphere','ball'} % a circle/ellipse in 2D
            radii = sqrt((1/(pi))*volume/prod(proportion)).*proportion;
            idx = (radii==0);
            radii(idx)=1;
            % First, specify a circular region.
            for d=1:2
               % region dimensions
               regionDims(d) = min(2*(ceil(min(radii)))+3,encompassingDims(d));
               % center coordinate of sphere template
               circleCenter(d) = (regionDims(d)+1)/2;
            end
            img = zeros(regionDims);
            % specify a exact circle with aspect ratio 1:1 with smallest radius
            disp('o Specifying circular template')
            circleRadius = min(radii) + 1;                      
            for i = 1:regionDims(1)
                for j = 1:regionDims(2)
                    if (i-circleCenter(1)+residualShift(1))^2 + (j-circleCenter(2)+residualShift(2))^2 <= circleRadius^2
                        img(i,j) = 1;
                    end
                end
            end
            disp('o Performing affine transformation.')            
            ratioMinDim = find(radii == min(radii),1,'first');
            Scale = [proportion(1)/(proportion(ratioMinDim)),0,0;0,proportion(2)/(proportion(ratioMinDim)),0;0,0,1];
            M = Scale'*Rotation';
            tshape = affine2d(M);
            RI = imref2d(size(img),1,1);
            newImg = imwarp(uint8(255*img),RI,tshape,method);
            newDims = size(newImg);
            padDims = ceil(0.5*(dimensions-newDims));         
            padDims(padDims<0) = 0;
            newImg = padarray(newImg,padDims);
            newDims = size(newImg);
            newDims = 2.*(ceil(0.5.*newDims));
            disp('o Building activation map.')            
            newImg = circshift(newImg,shiftVector-newDims/2);
            center = origin + center;
            for i = 1:dimensions(1)
                for j = 1:dimensions(2)
                    if newImg(i,j) > 0
                        if falloff ~= 0
                            activation_map(i,j) = (double(newImg(i,j))/255.0)*((1-minimum)*exp(-(((i-center(1))^2+(j-center(2))^2))*falloff)+minimum);
                        else
                            activation_map(i,j) = double(newImg(i,j))/255.0;
                        end
                    end
                end
            end
            
        case {'supercircle','superellipse','superellipsoid','superball','astroid','squircle','diamond'} % a superellipse in 2D
            radii = sqrt((gamma(1+2/n))/(4*(gamma(1+1/n))^2)*volume/prod(proportion)).*proportion;
            idx = (radii==0);
            radii(idx)=1;
            % First, specify a circular region.
            for d=1:2
               % region dimensions
               regionDims(d) = min(2*(ceil(min(radii)))+3,encompassingDims(d));
               % center coordinate of sphere template
               circleCenter(d) = (regionDims(d)+1)/2;
            end
            img = zeros(regionDims);
            % specify a exact circle with aspect ratio 1:1 with smallest radius
            disp('o Specifying circular template')
            circleRadius = min(radii) + 1;                      
            for i = 1:regionDims(1)
                for j = 1:regionDims(2)
                    if abs(i-circleCenter(1)+residualShift(1))^n + abs(j-circleCenter(2)+residualShift(2))^n <= circleRadius^n
                        img(i,j) = 1;
                    end
                end
            end
            disp('o Performing affine transformation.')            
            ratioMinDim = find(radii == min(radii),1,'first');
            Scale = [proportion(1)/(proportion(ratioMinDim)),0,0;0,proportion(2)/(proportion(ratioMinDim)),0;0,0,1];
            M = Scale'*Rotation';
            tshape = affine2d(M);
            RI = imref2d(size(img),1,1);
            newImg = imwarp(uint8(255*img),RI,tshape,method);
            newDims = size(newImg);
            padDims = ceil(0.5*(dimensions-newDims));         
            padDims(padDims<0) = 0;
            newImg = padarray(newImg,padDims);
            newDims = size(newImg);
            newDims = 2.*(ceil(0.5.*newDims));
            disp('o Building activation map.')            
            newImg = circshift(newImg,shiftVector-newDims/2);
            center = origin + center;
            for i = 1:dimensions(1)
                for j = 1:dimensions(2)
                    if newImg(i,j)> 0
                        if falloff ~= 0
                            activation_map(i,j) = (double(newImg(i,j))/255.0)*((1-minimum)*exp(-(((i-center(1))^2+(j-center(2))^2))*falloff)+minimum);
                        else
                            activation_map(i,j) = double(newImg(i,j))/255.0;
                        end
                    end
                end
            end
            
        case {'mask','manual'}                 
            % center is set of mask coordinates relative to origin
             center(:,1) = center(:,1) + origin(1);
             center(:,2) = center(:,2) + origin(2); 
             disp('o Building raw activation map.')              
             for i = 1:size(center,1)
                 x = center(i,1);
                 y = center(i,2);
                 activation_map(x,y) = 1;
             end      
             if falloff ~= 0
                maskCenter(1) = mean(center(:,1));
                maskCenter(2) = mean(center(:,2));
                maskRadius(1) = 2*std(center(:,1));
                maskRadius(2) = 2*std(center(:,2));
                proportion = [maskRadius(1) maskRadius(2)]/min(maskRadius);
                % define effect activation ellipse
                act.region = 'mask';
                act.volume = (4*pi)*prod(maskRadius);        
                act.center = maskCenter;  % left side of activation
                act.rotation = [0,0,0]; % degrees
                act.shape = 'ellipse';
                act.ratio = proportion;      % aspect ration
                act.falloff = falloff;  % parameterizes exponential falloff about center, in [0,1] 
                act.minimum = minimum;  % parameterizes exponential falloff minimum value in [0,1] 
                disp('o Building fading mask.') 
                falloffMask = STANCE_make_activation_map(dimensions, origin, act);
                activation_map = falloffMask.*activation_map;
             end              
             
        otherwise
            error('Unsupported shape ... check dimension?');
    end 
else
    switch shape
        
        case {'cube','box','cuboid','prism'}
            halfsides = 0.5*nthroot(volume/prod(proportion),3).*proportion;
            idx = (halfsides==0);
            halfsides(idx)=1;
            % First, specify the cube region 
            for d=1:3
               % region dimensions
               regionDims(d) = min(2*(ceil(min(halfsides)))+3,encompassingDims(d));
               % center coordinate of sphere template
               cubeCenter(d) = (regionDims(d)+1)/2;
            end
            img = zeros(regionDims);             
            % specify a exact cube with ratio 1:1:1 with smallest side
            disp('o Specifying cubical template')
            cubeHalfSide = min(halfsides) + 1;
            cubeCenter = (regionDims+1)/2;
            for i = 1:regionDims
                for j = 1:regionDims
                    for k = 1:regionDims
                        if abs(i-cubeCenter(1)+residualShift(1))<=cubeHalfSide && abs(j-cubeCenter(2)+residualShift(2))<=cubeHalfSide && abs(k-cubeCenter(3)+residualShift(3))<=cubeHalfSide
                            img(i,j,k) = 1;
                        end
                    end
                end
            end
            disp('o Performing affine transformation.')
            Scale = [proportion(1)/proportion(ratioMinDim),0,0,0;0,proportion(2)/proportion(ratioMinDim),0,0;0,0,proportion(3)/proportion(ratioMinDim),0;0,0,0,1];
            M = Scale'*Rx'*Ry'*Rz';
            tshape = affine3d(M);
            RI = imref3d(size(img),1,1,1);
            newImg = imwarp(uint8(255*img),RI,tshape,method);
            newDims = size(newImg);
            padDims = ceil(0.5*(dimensions-newDims));         
            padDims(padDims<0) = 0;
            newImg = padarray(newImg,padDims);
            newDims = size(newImg);
            newDims = 2.*(ceil(0.5.*newDims));
            disp('o Building activation map.')
            newImg = circshift(newImg,shiftVector-newDims/2);
            center = origin + center;            
            for i = 1:dimensions(1)
                for j = 1:dimensions(2)
                    for k = 1:dimensions(3)
                        if newImg(i,j,k) > 0
                            if falloff ~= 0
                                activation_map(i,j,k) = (double(newImg(i,j,k))/255.0)*((1-minimum)*exp(-((i-center(1))^2+(j-center(2))^2+(k-center(3))^2)*falloff)+minimum);
                            else
                                activation_map(i,j,k) = double(newImg(i,j,k))/255.0;
                            end
                        end
                    end
                end
            end
            
        case {'sphere','ellipsoid','ball'}
            radii = nthroot((3/(4*pi))*volume/prod(proportion),3).*proportion;
            idx = radii==0;
            radii(idx)=1;
            % firstly, set a sphere region to be specified. rl: region length
            for d=1:3
               % region dimensions
               regionDims(d) = min(2*(ceil(min(radii)))+3,encompassingDims(d));
               % center coordinate of sphere template
               sphereCenter(d) = (regionDims(d)+1)/2;
            end
            img = zeros(regionDims);            
            % specify a exact sphere with ratio 1:1:1 with smallest radius
            disp('o Specifying spherical template')
            % sphere radius
            sphereRadius = min(radii) + 1;  % add one to count center point
            for i = 1:regionDims(1)
                for j = 1:regionDims(2)
                    for k = 1:regionDims(3)
                        if (i-sphereCenter(1)+residualShift(1))^2 + (j-sphereCenter(2)+residualShift(2))^2 + (k-sphereCenter(3)+residualShift(3))^2 <= sphereRadius^2
                            img(i,j,k) = 1;
                        end
                    end
                end
            end          
            disp('o Performing affine transformation (this may take a while).')
            Scale = [proportion(1)/proportion(ratioMinDim),0,0,0;0,proportion(2)/proportion(ratioMinDim),0,0;0,0,proportion(3)/proportion(ratioMinDim),0;0,0,0,1];
            M = Scale'*Rx'*Ry'*Rz';
            tshape = affine3d(M);
            RI = imref3d(size(img),1,1,1);
            newImg = imwarp(uint8(255*img),RI,tshape,method);
            newDims = size(newImg);
            padDims = ceil(0.5*(dimensions-newDims));         
            padDims(padDims<0) = 0;
            newImg = padarray(newImg,padDims);
            newDims = size(newImg);
            newDims = 2.*(ceil(0.5.*newDims));
            disp('o Building activation map.')            
            newImg = circshift(newImg,shiftVector-newDims/2);            
            center = origin + center;
            for i = 1:dimensions(1)
                for j = 1:dimensions(2)
                    for k = 1:dimensions(3)
                        if newImg(i,j,k) > 0
                            if falloff ~= 0
                                activation_map(i,j,k) = (double(newImg(i,j,k))/255.0)*((1-minimum)*exp(-((i-center(1))^2+(j-center(2))^2+(k-center(3))^2)*falloff)+minimum);
                            else
                                activation_map(i,j,k) = (double(newImg(i,j,k))/255.0);
                            end
                        end
                    end
                end
            end                     

        case {'supersphere','superellipsoid','superball','superegg','astroid','squircle','diamond'}
            radii = nthroot(((3*r*n)/(8*beta(1/r,1/r)*beta(2/n,1/n)))*volume/prod(proportion),3).*proportion;
            idx = radii==0;
            radii(idx)=1;
            % firstly, set a sphere region to be specified. rl: region length
            for d=1:3
               % region dimensions
               regionDims(d) = min(2*(ceil(min(radii)))+3,encompassingDims(d));
               % center coordinate of sphere template
               sphereCenter(d) = (regionDims(d)+1)/2;
            end
            img = zeros(regionDims);            
            % specify a exact sphere with ratio 1:1:1 with smallest radius
            disp('o Specifying superspherical template')
            % sphere radius
            sphereRadius = min(radii) + 1;
            for i = 1:regionDims(1)
                for j = 1:regionDims(2)
                    for k = 1:regionDims(3)
                        if (abs(i-sphereCenter(1)+residualShift(1))^r + abs(j-sphereCenter(2)+residualShift(2))^r)^(n/r) + abs(k-sphereCenter(3)+residualShift(3))^n <= sphereRadius^n
                            img(i,j,k) = 1;
                        end
                    end
                end
            end          
            disp('o Performing affine transformation (this may take a while).')
            Scale = [proportion(1)/proportion(ratioMinDim),0,0,0;0,proportion(2)/proportion(ratioMinDim),0,0;0,0,proportion(3)/proportion(ratioMinDim),0;0,0,0,1];
            M = Scale'*Rx'*Ry'*Rz';
            tshape = affine3d(M);
            RI = imref3d(size(img),1,1,1);
            newImg = imwarp(uint8(255*img),RI,tshape,method);
            newDims = size(newImg);
            padDims = ceil(0.5*(dimensions-newDims));         
            padDims(padDims<0) = 0;
            newImg = padarray(newImg,padDims);
            newDims = size(newImg);
            newDims = 2.*(ceil(0.5.*newDims));
            disp('o Building activation map.')            
            newImg = circshift(newImg,shiftVector-newDims./2);            
            center = origin + center;
            for i = 1:dimensions(1)
                for j = 1:dimensions(2)
                    for k = 1:dimensions(3)
                        if newImg(i,j,k) > 0
                            if falloff ~= 0
                                activation_map(i,j,k) = (double(newImg(i,j,k))/255.0)*((1-minimum)*exp(-((i-center(1))^2+(j-center(2))^2+(k-center(3))^2)*falloff)+minimum);
                            else
                                activation_map(i,j,k) = (double(newImg(i,j,k))/255.0);
                            end
                        end
                    end
                end
            end      
            
        case {'mask','manual'}
            if falloff ~= 0
                warning('Activation fading for mask coordinates is not implemented yet')
            end
            % center is set of mask coordinates            
            center(:,1) = center(:,1) + origin(1);
            center(:,2) = center(:,2) + origin(2);
            center(:,3) = center(:,3) + origin(3);                                        
            for i = 1:volume
                x = center(i,1);
                y = center(i,2);
                k = center(i,3);
                activation_map(x,y,k) = 1;
            end
             if falloff ~= 0
                maskCenter(1) = mean(center(:,1));
                maskCenter(2) = mean(center(:,2));
                maskCenter(3) = mean(center(:,2));                
                maskRadius(1) = 2*std(center(:,1));
                maskRadius(2) = 2*std(center(:,2));
                maskRadius(3) = 2*std(center(:,3));                
                proportion = [maskRadius(1) maskRadius(2) maskRadius(3)]/min(maskRadius);
                % define effect activation ellipse
                act.region = 'mask';
                act.volume = (4*pi/3.0)*prod(maskRadius);        
                act.center = maskCenter;  % left side of activation
                act.rotation = [0,0,0]; % degrees
                act.shape = 'ellipse';
                act.ratio = proportion;      % aspect ration
                act.falloff = falloff;  % parameterizes exponential falloff about center, in [0,1] 
                act.minimum = minimum;  % parameterizes exponential falloff minimum value in [0,1] 
                disp('o Building fading mask.') 
                falloffMask = STANCE_make_activation_map(dimensions, origin, act);
                activation_map = falloffMask.*activation_map;
             end                          
             
        case {'atlas','parcellation'}
            disp('o Loading atlas ROI mask.') 

            if ~exist('STANCE.mat','file')
                if ~exist('..\STANCE.mat','file')
                   load('..\..\STANCE.mat'); 
                else
                   load('..\STANCE.mat');
                end
            else
                load('STANCE.mat');
            end

            filebase = region;
            if isempty(proportion)
               threshold = 0;
            else
               threshold = abs(proportion);
            end

            ROI_label = abs(volume);
            signFlag = sign(volume);
            switch filebase
                case 'aal'
                    filename = [STANCEroot,'/MNI/aal.nii'];
                    if ~exist(filename,'file')
                        filename = [filename,'.gz'];
                    end                    
                    [~,Y_atlas] = STANCE_load_volume(filename);
                    activation_map = (Y_atlas==ROI_label);
                    
                case 'brodmann'
                    filename = [STANCEroot,'/MNI/brodmann.nii'];
                    if ~exist(filename,'file')
                        filename = [filename,'.gz'];
                    end                    
                    [~,Y_atlas] = STANCE_load_volume(filename);
                    activation_map = (Y_atlas==ROI_label);      
                    
                case 'HarvardOxford'
                    if signFlag > 0 % cortical labels
                          if threshold == 0
                              filename = [STANCEroot,'/MNI/HarvardOxford-cort-maxprob-thr0-1mm.nii'];
                              if ~exist(filename,'file')
                                 filename = [filename,'.gz'];
                              end     
                          elseif threshold == 25 || threshold == 0.25
                              filename = [STANCEroot,'/MNI/HarvardOxford-cort-maxprob-thr25-1mm.nii'];
                              if ~exist(filename,'file')
                                 filename = [filename,'.gz'];
                              end
                          elseif threshold == 50 || threshold == 0.5
                              filename = [STANCEroot,'/MNI/HarvardOxford-cort-maxprob-thr50-1mm.nii'];                     
                              if ~exist(filename,'file')
                                 filename = [filename,'.gz'];
                              end
                          else
                              error('Unsupported probability threshold')
                          end
                    else % sub-cortical labels
                          if threshold == 0
                              filename = [STANCEroot,'/MNI/HarvardOxford-sub-maxprob-thr0-1mm.nii.gz'];
                          elseif threshold == 25 || threshold == 0.25
                              filename = [STANCEroot,'/MNI/HarvardOxford-sub-maxprob-thr25-1mm.nii.gz'];
                          elseif threshold == 50 || threshold == 0.5
                              filename = [STANCEroot,'/MNI/HarvardOxford-sub-maxprob-thr50-1mm.nii.gz'];                     
                          else
                              error('Unsupported probability threshold')
                          end                          
                    end
                    [~,Y_atlas] = STANCE_load_volume(filename);
                    Y_atlas = Y_atlas(2:end,2:end,2:end);
                    activation_map = (Y_atlas==ROI_label); 
                    
                case 'Craddock'
                    if threshold == 200
                        filename = [STANCEroot,'/MNI/CC200ROI_tcorr05_2level_all.nii'];
                        if ~exist(filename,'file')
                            filename = [filename,'.gz'];
                        end
                    elseif threshold == 400
                        filename = [STANCEroot,'/MNI/CC400ROI_tcorr05_2level_all.nii'];
                        if ~exist(filename,'file')
                            filename = [filename,'.gz'];
                        end                        
                    else
                        error('Unsupported number of parcellations')
                    end
                    [~,Y_temp] = STANCE_load_volume(filename);
                    Y_atlas = zeros(188,224,184);
                    for i = 1:size(Y_temp,1)
                        for j = 1:size(Y_temp,2)
                            for k = 1:size(Y_temp,3)
                                Y_atlas(4*i-3:4*i,4*j-3:4*j,4*k-3:4*k) = Y_temp(i,j,k);
                            end
                        end
                    end
                    Y_atlas = Y_atlas(186:-1:6,7:223,4:184);
                    activation_map = (Y_atlas==ROI_label); 
                    
                case 'Brainnetome'
                    if threshold == 25
                        filename = [STANCEroot,'/MNI/BrainnetomeAtlas_BNA_MPM_thr25_1.25mm.nii'];
                        if ~exist(filename,'file')
                            filename = [filename,'.gz'];
                        end                      
                    else
                        error('Unsupported probability threshold')
                    end
                    [~,Y_temp] = STANCE_load_volume(filename);
                    Y_temp(:,174,:) = 0;
                    Y_atlas = zeros(181,218,183);
                    for i = 1:181
                        for j = 1:217
                            for k = 1:181
                                Y_atlas(i,j,k) = Y_temp(round(0.8*i),round(0.8*j),round(0.8*k));
                            end
                        end
                    end
                    Y_atlas = Y_atlas(181:-1:1,2:218,3:183);
                    activation_map = (Y_atlas==ROI_label); 
                    
                otherwise
                       filename = filebase;
                       [~,Y_atlas] = STANCE_load_volume(filename);
                       if sum(size(Y_atlas) > 0) > 3
                           activation_map = Y_atlas(:,:,:,ROI_label);
                       else
                           activation_map = (Y_atlas==ROI_label);
                       end
            end
            if falloff ~= 0
                maskCenter = center_of_mass(activation_map);               
                maskSide(1) = sum(sum(sum(activation_map,2),3)>0);
                maskSide(2) = sum(sum(sum(activation_map),3)>0);
                maskSide(3) = sum(sum(sum(activation_map),2)>0);               
                proportion = [maskSide(1) maskSide(2) maskSide(3)]/min(maskSide);
                % define effect activation ellipse
                act.region = 'mask';
                act.volume = prod(maskSide);        
                act.center = maskCenter';  % left side of activation
                act.rotation = [0,0,0]; % degrees
                act.shape = 'box';
                act.ratio = proportion;      % aspect ration
                act.falloff = falloff;  % parameterizes exponential falloff about center, in [0,1] 
                act.minimum = minimum;  % parameterizes exponential falloff minimum value in [0,1] 
                disp('o Building fading mask.') 
                falloffMask = STANCE_make_activation_map(dimensions, origin, act);
                activation_map = falloffMask.*activation_map;
            end                          

        case {'data'}
            disp('o Loading data from file.')             

            if exist('../../STANCE.mat','file')
                load('../../STANCE.mat')
            elseif exist('../STANCE.mat','file')
                load('../STANCE.mat')
            elseif exist('STANCE.mat','file')
                    load('STANCE.mat')
            end            
            
            [V_MNI,~] = STANCE_load_volume(filenameMNI);
            
            activation_map = STANCE_load_map(region,V_MNI,proportion);
            
        otherwise
            error('Unsupported shape ... check dimension?');            
    end
end

activation_map = double(activation_map);

end