
% BFILTER3 Three dimensional bilateral filtering.
%    This function implements 3-D bilateral filtering using
%    the method outlined in:
%
%       C. Tomasi and R. Manduchi. Bilateral Filtering for 
%       Gray and Color Images. In Proceedings of the IEEE 
%       International Conference on Computer Vision, 1998. 
%
%    B = bfilter3(A,W,SIGMA) performs 3-D bilateral filtering
%    for the grayscale or color image A. A should be a double
%    precision matrix of size LxNxMx1 or LxNxMx3 (i.e., grayscale
%    or color images, respectively) with normalized values in
%    the closed interval [0,1]. The half-size of the Gaussian
%    bilateral filter window is defined by W. The standard
%    deviations of the bilateral filter are given by SIGMA,
%    where the spatial-domain standard deviation is given by
%    SIGMA(1) and the intensity-domain standard deviation is
%    given by SIGMA(2).
%
% Author: adapted for use in 3D Volumes by Kevin Matlock.
% Original 2D Filter By: Douglas R. Lanman, Brown University, 
% September 2006. dlanman@brown.edu, http://mesh.brown.edu/dlanman
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-process input and select appropriate filter.
function B = bfilter3(A,w,sigma)

% Verify that the input image exists and is valid.
if ~exist('A','var') || isempty(A)
   error('Input image A is undefined or invalid.');
end
if ~isfloat(A) || ~sum([1,3] == size(A,4)) || ...
      min(A(:)) < 0 || max(A(:)) > 1
   error(['Input image A must be a double precision ',...
          'matrix of size LxNxMx1 or LxNxMx3 on the closed ',...
          'interval [0,1].']);      
end

% Verify bilateral filter window size.
if ~exist('w','var') || isempty(w) || ...
      numel(w) ~= 1 || w < 1
   w = 5;
end
w = ceil(abs(w));

% Verify bilateral filter standard deviations.
if ~exist('sigma','var') || isempty(sigma) || numel(sigma) ~= 2
   sigma = [3 0.1];
else
   sigma = [abs(sigma(1)) abs(sigma(2))];
end

% Apply either grayscale or color bilateral filtering.
if size(A,4) == 1
   B = bfltGray3(A,w,sigma(1),sigma(2));
else
   B = bfltColor3(A,w,sigma(1),sigma(2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements bilateral filtering for grayscale images.
function B = bfltGray3(A,w,sigma_d,sigma_r)

% Pre-compute Gaussian domain weights.
[X,Y,Z] = meshgrid(-w:w,-w:w, -w:w);
G = exp(-(X.^2+Y.^2+Z.^2)/(2*sigma_d^2));

% Create waitbar.
h = waitbar(0,'Applying bilateral filter...');
set(h,'Name','Bilateral Filter Progress');

% Apply bilateral filter.
dim = size(A);
B = zeros(dim);
for i = 1:dim(1)
   for j = 1:dim(2)
       for k = 1:dim(3)
      
         % Extract local region.
         iMin = max(i-w,1);
         iMax = min(i+w,dim(1));
         jMin = max(j-w,1);
         jMax = min(j+w,dim(2));
         kMin = max(k-w,1);
         kMax = min(k+w,dim(3));
         I = A(iMin:iMax,jMin:jMax,kMin:kMax);
      
         % Compute Gaussian range weights.
         H = exp(-(I-A(i,j,k)).^2/(2*sigma_r^2));
         
         % Calculate bilateral filter response.
         F = H.*G((iMin:iMax)-i+w+1,(jMin:jMax)-j+w+1,(kMin:kMax)-k+w+1);
         norm_F = sum(F(:));
         B(i,j,k) = sum(F(:).*I(:))/norm_F;
    
       end
   end
   waitbar(i/dim(1));
end

% Close waitbar.
close(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements bilateral filter for color images.
function B = bfltColor3(A,w,sigma_d,sigma_r)

A = colorspace3d(A, 1);

% Pre-compute Gaussian domain weights.
[X,Y,Z] = meshgrid(-w:w,-w:w, -w:w);
G = exp(-(X.^2+Y.^2+Z.^2)/(2*sigma_d^2));

% Rescale range variance (using maximum luminance).
sigma_r = 100*sigma_r;

% Create waitbar.
h = waitbar(0,'Applying bilateral filter...');
set(h,'Name','Bilateral Filter Progress');

% Apply bilateral filter.
dim = size(A);
B = zeros(dim);
for i = 1:dim(1)
   for j = 1:dim(2)
       for k = 1:dim(3)
      
         % Extract local region.
         iMin = max(i-w,1);
         iMax = min(i+w,dim(1));
         jMin = max(j-w,1);
         jMax = min(j+w,dim(2));
         kMin = max(k-w,1);
         kMax = min(k+w,dim(3));
         I = A(iMin:iMax,jMin:jMax,kMin:kMax,:);
      
         % Compute Gaussian range weights.
         dL = I(:,:,:,1)-A(i,j,k,1);
         da = I(:,:,:,2)-A(i,j,k,2);
         db = I(:,:,:,3)-A(i,j,k,3);
         H = exp(-(dL.^2+da.^2+db.^2)/(2*sigma_r^2));
      
         % Calculate bilateral filter response.
         F = H.*G((iMin:iMax)-i+w+1,(jMin:jMax)-j+w+1,(kMin:kMax)-k+w+1);
         norm_F = sum(F(:));
         B(i,j,k,1) = sum(sum(sum(F.*I(:,:,:,1))))/norm_F;
         B(i,j,k,2) = sum(sum(sum(F.*I(:,:,:,2))))/norm_F;
         B(i,j,k,3) = sum(sum(sum(F.*I(:,:,:,3))))/norm_F;
                
       end
   end
   waitbar(i/dim(1));
end

B = colorspace3d(B, 0);

% Close waitbar.
close(h);