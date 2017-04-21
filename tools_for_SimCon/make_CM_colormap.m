function map = make_CM_colormap(Ncolors,zerocolor)
% Generates correlation matrix color map: 
% [dark blue -> blue -> cyan -> white -> yellow -> red -> dark red]
%
% Author: Jason E. Hill, Ph.D.
% Date: 31 AUG 2016

% check that the integer N has 8 as a factor

if nargin < 2
    zerocolor = [1; 1; 1];
else
    if ischar(zerocolor)
        zerocolor = squeeze(double(label2rgb(0, jet, zerocolor)))/255;
    end
end
if mod(Ncolors,8) == 0
    % good to go
else
    warning('Number of color channels must be a positive multiple of 8! Rounding ...');
    Ncolors = round(abs(Ncolors)/8)*8;
    if Ncolors < 16
        Ncolors = 16;
    end
end

% initialize 
map = zeros(Ncolors,3);
map(((Ncolors+8)/8):(Ncolors/4)  ,3) = 1;
map(((Ncolors+4)/4):(Ncolors/4)  ,2) = 1;
map(((3*Ncolors+4)/4):(7*Ncolors/8),1) = 1;
% dark blue to blue transition
map(1:(Ncolors/8)            ,3) = 0.5+0.5*(0:(8/(Ncolors-8)):1);
% blue to cyan transition
map(((Ncolors+8)/8):(Ncolors/4)    ,2) = (0:(8/(Ncolors-8)):1);
% cyan to zerocolor transition
map(((Ncolors+4)/4):(Ncolors/2)    ,1) = zerocolor(1)*(0:(4/(Ncolors-4)):1);
map(((Ncolors+4)/4):(Ncolors/2)    ,2) = 1-(1-zerocolor(2))*(0:(4/(Ncolors-4)):1);
map(((Ncolors+4)/4):(Ncolors/2)    ,3) = 1-(1-zerocolor(3))*(0:(4/(Ncolors-4)):1);
% zerocolor to yellow transition
map(((Ncolors+2)/2):(3*Ncolors/4)  ,1) = 1-(1-zerocolor(1))*(1:(-4/(Ncolors-4)):0);
map(((Ncolors+2)/2):(3*Ncolors/4)  ,2) = 1-(1-zerocolor(2))*(1:(-4/(Ncolors-4)):0);
map(((Ncolors+2)/2):(3*Ncolors/4)  ,3) = zerocolor(3)*(1:(-4/(Ncolors-4)):0);
% yellow to red transition
map(((3*Ncolors+4)/4):(7*Ncolors/8),2) = (1:(-8/(Ncolors-8)):0);
% red dark red transition
map(((7*Ncolors+8)/8):Ncolors      ,1) = (1:(-8/(Ncolors-8)):0);

end