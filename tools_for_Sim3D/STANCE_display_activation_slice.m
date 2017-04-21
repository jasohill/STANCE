function h = STANCE_display_activation_slice(Y_brain,Y_activation,slice,direction,origin,GUIflag)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% Overlays activation on brain for slice in dimension.
%
% Jason E. Hill
% STANCE_display_activation_slice.m      updated     28 SEP 2016

if nargin < 6
    GUIflag = false;
end
if nargin < 5
    origin = [0 0 0];
end
if nargin < 4
    direction = 3;  % axial view is default
end
if nargin < 3       % mid-point is default
    slice = round(size(Y_brain,3)/2);
end

if isempty(direction) % show all three
    if GUIflag
        h3 = STANCE_display_activation_slice(Y_brain,Y_activation,slice,3,origin,GUIflag);
        h1 = STANCE_display_activation_slice(Y_brain,Y_activation,slice,1,origin,GUIflag);
        h2 = STANCE_display_activation_slice(Y_brain,Y_activation,slice,2,origin,GUIflag);
    else
        GUIflag = true; % use to turn off figure
        [X,Y,Z] = size(Y_brain);
        height = Y+Z+70;
        width  = X+Y+40;
        positionVectorC = [20/width,    20/height, X/width, Z/height];
        positionVectorS = [(X+20)/width,20/height, Y/width, Z/height];        
        positionVectorA = [20/width,(Z+20)/height, X/width, Y/height];
        h = figure;
        subplot('Position',positionVectorC)
        STANCE_display_activation_slice(Y_brain,Y_activation,slice,2,origin,GUIflag);       
        subplot('Position',positionVectorS)
        STANCE_display_activation_slice(Y_brain,Y_activation,slice,1,origin,GUIflag);
        subplot('Position',positionVectorA)
        STANCE_display_activation_slice(Y_brain,Y_activation,slice,3,origin,GUIflag);        
    end  
elseif direction <0.5 || direction>3.5
    warning('Unrecognized direction, showing all three.')
    h3 = STANCE_display_activation_slice(Y_brain,Y_activation,slice,3,origin,GUIflag);
    h1 = STANCE_display_activation_slice(Y_brain,Y_activation,slice,1,origin,GUIflag);
    h2 = STANCE_display_activation_slice(Y_brain,Y_activation,slice,2,origin,GUIflag);   
else 
    direction = round(abs(direction));
    if isempty(slice)     % choose slice with maximum activations
        if direction == 1 % supress mid-sagittal activations
           Y_activation(floor(0.495*size(Y_activation,1)):ceil(0.505*size(Y_activation,1)),:,:) = 0;
        end
        if sum(Y_activation(:)) == 0
            slice = round(size(Y_activation,direction)/2);
        else
            if direction == 1     % sagittal
                [~,J_max] = max( sum(sum(Y_activation,2),3));
                slice = J_max(1);   
            elseif direction == 2 % coronal
                [~,J_max] = max(sum(sum(Y_activation),3));
                slice = J_max(1);
            else                  % axial
                [~,J_max] = max(sum(sum(Y_activation)));
                slice = J_max(1);                
            end   
        end
    end    
    
if isempty(Y_brain)

if direction == 1     % sagittal 
    J = imrotate(squeeze(uint8(255*Y_activation(slice,:,:))),90);
elseif direction == 2 % coronal
    J = imrotate(squeeze(uint8(255*Y_activation(:,slice,:))),90);    
else                  % axial
    J = imrotate(uint8(255*Y_activation(:,:,slice)),90);
end    

if size(J,2) == size(J,1)
    C5  = round(3*size(J,2)/181.0 +2);
    C8  = round(6*size(J,1)/181.0 +2);
    C10 = round(8*size(J,1)/181.0 +2);
    C16 = round(14*size(J,2)/181.0+2);
    C25 = round(23*size(J,1)/181.0+2);
elseif size(J,2) > size(J,1)
    C5  = round(3*size(J,2)/217.0 +2);
    C8  = round(6*size(J,1)/181.0 +2);
    C10 = round(8*size(J,1)/181.0 +2);
    C16 = round(14*size(J,2)/217.0+2);
    C25 = round(23*size(J,1)/181.0+2);
else
    C5  = round(3*size(J,2)/181.0 +2);
    C8  = round(6*size(J,1)/217.0 +2);
    C10 = round(8*size(J,1)/217.0 +2);
    C16 = round(14*size(J,2)/181.0+2);
    C25 = round(23*size(J,1)/217.0+2);   
end

if GUIflag
h = imshow(J,[]);
if direction == 1
    text(C5,C10,num2str(round(slice-origin(1))),'Color','red','FontSize',14)     
    text(C5,C25,'P','Color','green','FontSize',14)
    text(size(J,2)-C16,C25,'A','Color','green','FontSize',14) 
    text(size(J,2)-C16,C10,'S','Color','blue','FontSize',14)
    text(size(J,2)-C16,size(J,1)-C8,'I','Color','blue','FontSize',14) 
elseif direction == 2
    text(C5,C10,num2str(round(slice-origin(2))),'Color','green','FontSize',14)         
    text(size(J,2)-C16,C10,'S','Color','blue','FontSize',14)
    text(size(J,2)-C16,size(J,1)-C8,'I','Color','blue','FontSize',14) 
    text(C5,C25,'L','Color','red','FontSize',14)
    text(size(J,2)-C16,C25,'R','Color','red','FontSize',14)
else
    text(C5,C10,num2str(round(slice-origin(3))),'Color','blue','FontSize',14)    
    text(C5,C25,'L','Color','red','FontSize',14)
    text(size(J,2)-C16,C25,'R','Color','red','FontSize',14)
    text(size(J,2)-C16,C10,'A','Color','green','FontSize',14)
    text(size(J,2)-C16,size(J,1)-C8,'P','Color','green','FontSize',14)     
end       
else
figure,    
h = imshow(J,[]);
if direction == 1
    text(C5,C10,num2str(round(slice-origin(1))),'Color','red','FontSize',14)     
    text(C5,C25,'P','Color','green','FontSize',14)
    text(size(J,2)-C16,C25,'A','Color','green','FontSize',14) 
    text(size(J,2)-C16,C10,'S','Color','blue','FontSize',14)
    text(size(J,2)-C16,size(J,1)-C8,'I','Color','blue','FontSize',14) 
elseif direction == 2
    text(C5,C10,num2str(round(slice-origin(2))),'Color','green','FontSize',14)         
    text(size(J,2)-C16,C10,'S','Color','blue','FontSize',14)
    text(size(J,2)-C16,size(J,1)-C8,'I','Color','blue','FontSize',14) 
    text(C5,C25,'L','Color','red','FontSize',14)
    text(size(J,2)-C16,C25,'R','Color','red','FontSize',14)
else
    text(C5,C10,num2str(round(slice-origin(3))),'Color','blue','FontSize',14)    
    text(C5,C25,'L','Color','red','FontSize',14)
    text(size(J,2)-C16,C25,'R','Color','red','FontSize',14)
    text(size(J,2)-C16,C10,'A','Color','green','FontSize',14)
    text(size(J,2)-C16,size(J,1)-C8,'P','Color','green','FontSize',14)     
end          
end
    
else
Ymax = max(Y_brain(:));

if direction == 1      % axial
    I = imrotate(squeeze(uint8(255*(Y_brain(slice,:,:)/Ymax))),90);
    J = imrotate(squeeze(uint8(255*Y_activation(slice,:,:))),90);
elseif direction == 2 % sagittal
    I = imrotate(squeeze(uint8(255*(Y_brain(:,slice,:)/Ymax))),90);
    J = imrotate(squeeze(uint8(255*Y_activation(:,slice,:))),90);    
else                  % coronal
    I = imrotate(uint8(255*(Y_brain(:,:,slice)/Ymax)),90);
    J = imrotate(uint8(255*Y_activation(:,:,slice)),90);
end

I(J>0) = 0;

I_rgb(:,:,1) = I + uint8(255*(J>0));
I_rgb(:,:,2) = I + J;
I_rgb(:,:,3) = I;

if size(J,2) == size(J,1)
    C5  = round(3*size(J,2)/181.0+2);
    C8  = round(4*size(J,1)/181.0+4);
    C10 = round(5*size(J,1)/181.0+5);
    C16 = round(8*size(J,2)/181.0+8);
    C25 = round(13*size(J,1)/181.0+12);
elseif size(J,2) > size(J,1)
    C5  = round(3*size(J,2)/217.0+2);
    C8  = round(4*size(J,1)/181.0+4);
    C10 = round(5*size(J,1)/181.0+5);
    C16 = round(8*size(J,2)/217.0+8);
    C25 = round(13*size(J,1)/181.0+12);
else
    C5  = round(3*size(J,2)/181.0+2);
    C8  = round(4*size(J,1)/217.0+4);
    C10 = round(5*size(J,1)/217.0+5);
    C16 = round(8*size(J,2)/181.0+8);
    C25 = round(13*size(J,1)/217.0+12);   
end
    
if GUIflag
h = imshow(I_rgb);
if direction == 1
    text(C5,C10,num2str(round(slice-origin(1))),'Color','red','FontSize',14)     
    text(C5,C25,'P','Color','green','FontSize',14)
    text(size(I,2)-C16,C25,'A','Color','green','FontSize',14) 
    text(size(I,2)-C16,C10,'S','Color','blue','FontSize',14)
    text(size(I,2)-C16,size(I,1)-C8,'I','Color','blue','FontSize',14) 
elseif direction == 2
    text(C5,C10,num2str(round(slice-origin(2))),'Color','green','FontSize',14)         
    text(size(I,2)-C16,C10,'S','Color','blue','FontSize',14)
    text(size(I,2)-C16,size(I,1)-C8,'I','Color','blue','FontSize',14) 
    text(C5,C25,'L','Color','red','FontSize',14)
    text(size(I,2)-C16,C25,'R','Color','red','FontSize',14)
else
    text(C5,C10,num2str(round(slice-origin(3))),'Color','blue','FontSize',14)    
    text(C5,C25,'L','Color','red','FontSize',14)
    text(size(I,2)-C16,C25,'R','Color','red','FontSize',14)
    text(size(I,2)-C16,C10,'A','Color','green','FontSize',14)
    text(size(I,2)-C16,size(I,1)-C8,'P','Color','green','FontSize',14)     
end           
else
figure,
h = imshow(I_rgb);
if direction == 1
    text(C5+1,C10+2,num2str(round(slice-origin(1))),'Color','red','FontSize',13)     
    text(C5+1,C25+5,'P','Color','green','FontSize',13)
    text(size(I,2)-C16-3,C25+5,'A','Color','green','FontSize',13) 
    text(size(I,2)-C16-3,C10+2,'S','Color','blue','FontSize',13)
    text(size(I,2)-C16-3,size(I,1)-C8-2,'I','Color','blue','FontSize',13) 
elseif direction == 2
    text(C5+1,C10+2,num2str(round(slice-origin(2))),'Color','green','FontSize',13)         
    text(size(I,2)-C16-3,C10,'S','Color','blue','FontSize',13)
    text(size(I,2)-C16-3,size(I,1)-C8-2,'I','Color','blue','FontSize',13) 
    text(C5+1,C25+5,'L','Color','red','FontSize',13)
    text(size(I,2)-C16-3,C25+5,'R','Color','red','FontSize',13)
else
    text(C5+1,C10+2,num2str(round(slice-origin(3))),'Color','blue','FontSize',13)    
    text(C5+1,C25+5,'L','Color','red','FontSize',13)
    text(size(I,2)-C16-3,C25+5,'R','Color','red','FontSize',13)
    text(size(I,2)-C16-3,C10+2,'A','Color','green','FontSize',13)
    text(size(I,2)-C16-3,size(I,1)-C8-2,'P','Color','green','FontSize',13)     
end       
end
end
end
end

