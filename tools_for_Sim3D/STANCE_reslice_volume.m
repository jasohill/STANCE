function [V_new,Y_new] = STANCE_reslice_volume(V_old,scan,sumThreshold)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% STANCE_reslice_volume generates a resliced volume V_new of V_old according to 
% the fMRI scan parameters that specify the functional space:
%   sumThreshold if an integer this quantity sets the requirement for a slice to be
%   included; if an array it specifies the slice limits
%
%        V = same format as SPM12 header
%        Y = MRI data array
%
%        scan.voxel.size: 3x1 list of voxel physical dimensions [mm]
%                        (default) = [3 3 3]
%        scan.voxel.dimensions: 3x1 list of the number of voxels per spatial direction
%                        (default) = [64 64 NaN]
%                        NOTE: dimensions(3) = NaN used indicate auto
%                              determine setting
%        scan.voxel.spacing: 3x1 list of inter-voxel spacing [mm]
%        scan.tiltAngle: tilt angle in YZ plane away from nose [degrees]
%                        (default) = 15
%
%        sumThreshold: (option 1) numerical threshold required upon a slice 
%                                 for inclusion in the resulting 3D volume.
%                      (option 2) a 2 element array where
%                                 sumThreshold(1) = the lower slice index
%                                 to use after the affine transformation
%                                 sumThreshold(1) = the upper slice index
%                                 to use after the affine transformation   
%
% author: Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
% STANCE_reslice_volume.m      updated     28 SEP 2016

%% Input argument handling and initialization

if nargin < 3
    sumThreshold   = 100;
end
if length(sumThreshold) == 2
    sliceLimitLower = sumThreshold(1);
    sliceLimitUpper = sumThreshold(2);
    sliceLimitFlag  = true;
else
    sliceLimitFlag  = false;
end    
if nargin < 2
    scan = [];
end

dimMax = 320; % limit any dimension to this for sake of memory consumption/error handling

%% Load scan info 

if (nargin < 2) || isempty(scan)
    voxelSize    = [3 3 3];
    new_dims     = [64 64 NaN];
    tiltAngle    = 15; % degrees
    voxelSpacing = [0 0 0.6];
else
    voxelSize    = scan.voxel.size;
    new_dims     = scan.voxel.matrix;
    voxelSpacing = scan.voxel.spacing;    
    tiltAngle    = scan.tiltAngle;
end
if isnan(new_dims(3))
    sliceNbrFlag = false;
else
    sliceNbr = new_dims(3);
    sliceNbrFlag = true;
end

dX = voxelSize(1) + voxelSpacing(1);
dY = voxelSize(2) + voxelSpacing(2);
dZ = voxelSize(3) + voxelSpacing(3);

%% Load header info

old_dims = [V_old.dim(1), V_old.dim(2), V_old.dim(3)]';
origin = [abs(V_old.mat(1,4)),abs(V_old.mat(2,4)),abs(V_old.mat(3,4))]';
dX_old = abs(V_old.mat(1,1));
dY_old = abs(V_old.mat(2,2));
dZ_old = abs(V_old.mat(3,3));

if sum(old_dims-origin) < double(sum(old_dims))/30.0
    error('The origin is not near the center of the image!')
end
    
%% determine necessary values for affine transformation

Xc = (old_dims(1) + 1)/2;
Yc = (old_dims(2) + 1)/2;
Zc = (old_dims(3) + 1)/2;

Ct = cosd(tiltAngle);
St = sind(tiltAngle);

nZ = ceil((Zc*Ct + Yc*St)/(dZ/dZ_old));
Zc = origin(3);

%% reslice slice by slice according to affine transformation

startFlag  = true;
stopFlag   = false;

nMax = (2*nZ+1);
nMax = min(nMax,dimMax);

hw = waitbar(0,'Reslicing volume...');
for n = 1:nMax
    if ~stopFlag
        % set rotation about x-axis thru point (x,Yc,Zc)
        % keeping in plane dimensions the same
        %  -> new re-sliced centers along the line: Zk = Zc + T(Yk - Yc),
        %  where the re-sliced center points are at
        %    Ykc = Yc + k*dZ*St
        %    Zkc = Zc + k*dZ*Ct
        %     k = -nZ ... 0 ... nZ
        %     n = k + (nZ + 1)
        k = n - (nZ+1);
        % re-sliced origin:
        Yk = Yc*(1 - Ct) - k*dZ*St;
        Zk = Zc - Yc*St + k*dZ*Ct; 
        Affine(:,:,n) = [[1 0   0   0]; ...
                         [0 Ct -St  Yk]; ...
                         [0 St  Ct  Zk]; ...                
                         [0 0   0   1]];           %#ok<*AGROW>
        affine_dims = [old_dims(1), old_dims(2)];
        AffineSlice = spm_slice_vol(V_old,Affine(:,:,n),affine_dims,1);    
        
        slicesSums(n) = squeeze(sum(sum(AffineSlice)));
%        figure, imshow(AffineSlice);
        if  ~sliceLimitFlag
            if sliceNbrFlag
                AffineSlices(:,:,n) = AffineSlice;
            elseif (slicesSums(n) > sumThreshold) && startFlag && ~stopFlag
                startFlag = false;
                AffineSlices(:,:,1) = AffineSlice;
                slice_n = 1;
                sliceLimitLower = n;
            elseif (slicesSums(n) > sumThreshold) && ~startFlag && ~stopFlag
                slice_n = slice_n + 1;
                AffineSlices(:,:,slice_n) = AffineSlice;
%                figure, imshow(sA);
            elseif (slicesSums(n) < sumThreshold) && ~startFlag && ~stopFlag
                stopFlag = true;
                sliceLimitUpper = n-1;
            end
        else
            if n == sliceLimitLower
                slice_n = 1;
                AffineSlices(:,:,1) = AffineSlice;
                startFlag = false;
            elseif (n > sliceLimitLower)&&(n < sliceLimitUpper)
                slice_n = slice_n + 1;
                AffineSlices(:,:,slice_n) = AffineSlice;
            end
            if n == sliceLimitUpper
                slice_n = slice_n + 1;
                AffineSlices(:,:,slice_n) = AffineSlice;
                stopFlag = true;
            end 
        end
    end
    waitbar(n/nMax,hw)
end
if (sliceNbrFlag && ~sliceLimitFlag)
    Nslices = length(slicesSums);
    % add small number to guide selection of empty slices toward the center
    % prefering the lower slices
    slicesSums(1:floor(0.5*Nslices)) = slicesSums(1:floor(0.5*Nslices)) + 0.01*(1:floor(0.5*Nslices));
    slicesSums(ceil(0.5*Nslices):end) = slicesSums(ceil(0.5*Nslices):end) + 0.01*(fliplr(ceil(0.5*Nslices):Nslices)-(ceil(0.5*Nslices)-1.5));
    [~, idxSort]  = sort(slicesSums,'descend');
    if sliceNbr < Nslices
        sliceLimit = min(idxSort(1:sliceNbr));
    else
        sliceLimit = 1;
    end
    sliceLimitLower = sliceLimit;
    if sliceLimit + sliceNbr - 1 > nMax          
        sliceLimitUpper = nMax;
        sliceLimitLower = sliceLimit - (nMax - (sliceLimit + sliceNbr - 1));
        if sliceLimitLower < 1
            sliceLimitLower = 1;
            sliceLimitUpper = sliceNbr;
            AffineSlices = padarray(AffineSlices,[0 0 (sliceNbr-nMax)],0,'post');
            sliceNbrWarning = sprintf('To accomodate the requested number of slices, zero-padding by %d slices!',(sliceNbr-nMax));
            warning(sliceNbrWarning); %#ok<*SPWRN>
        end
    else
        sliceLimitUpper = sliceLimit + sliceNbr - 1;
    end
    AffineSlices = AffineSlices(:,:,sliceLimitLower:sliceLimitUpper);
    slice_n = size(AffineSlices,3);
end

affine_dims = size(AffineSlices);
affine_dims(1) = min(affine_dims(1),dimMax);
affine_dims(2) = min(affine_dims(2),dimMax);

AffineOrigin(1) = origin(1);
AffineOrigin(2) = Yc + (origin(2)-Yc)*Ct + (origin(3)-Zc)*St;
AffineOrigin(3) = (nZ - sliceLimitLower + 2) + (origin(3) - Zc)*Ct/dZ + (Yc - origin(2))*St/dZ;

%% reslice to target scan dimensions

affineImageCenter(1) = 0.5*(affine_dims(1)+1);
affineImageCenter(2) = 0.5*(affine_dims(2)+1);

newImageCenter(1) = 0.5*(new_dims(1)+1);
newImageCenter(2) = 0.5*(new_dims(2)+1);

newOrigin(1) = newImageCenter(1) + (AffineOrigin(1) - Xc)/dX;
newOrigin(2) = newImageCenter(2) + (AffineOrigin(2) - Yc)/dY;
newOrigin(3) = AffineOrigin(3);

[YA,XA] = meshgrid(1:affine_dims(2),1:affine_dims(1));
[YN,XN] = meshgrid(affineImageCenter(2)-0.5*dY*(new_dims(2)-1):dY:affineImageCenter(2)+0.5*dY*(new_dims(2)-1),...
                   affineImageCenter(1)-0.5*dX*(new_dims(1)-1):dX:affineImageCenter(1)+0.5*dX*(new_dims(1)-1));
               
newData = zeros(new_dims(1),new_dims(2),slice_n);
%size(newData)
for n = 1:slice_n
    newData(:,:,n) = interp2(YA,XA,AffineSlices(:,:,n),YN,XN);
%    figure, imshow(newData(:,:,n));
end

%% save to new volume structure
fileName = V_old.fname;
fileExtension = fileName(end-2:end);
if strcmpi(fileExtension,'.gz')
    fileExtension = fileName(end-6:end-3);
    fileBase = fileName(1:end-7);
else
    fileBase = fileName(1:end-4);
    fileExtension = fileName(end-3:end);
end
dataFileNameNIIout = [fileBase,'_resliced',fileExtension];

V_new.fname = dataFileNameNIIout;
V_new.dim = size(newData);
V_new.dt = V_old.dt;
V_new.pinfo = V_old.pinfo;
V_new.n = [1,1];
V_new.descrip = V_old.descrip;
V_new.private = V_old.private;
Y_new = newData;
Y_new(isnan(Y_new)) = 0;

if length(sumThreshold) == 1
    V_new.sliceLimitLower = sliceLimitLower;
    V_new.sliceLimitUpper = sliceLimitUpper;
end

% generate new world coordinate matrix
sign_mat = sign(V_old.mat);
new_mat = zeros(4);
new_mat(1,1) = sign_mat(1,1)*dX; 
new_mat(2,2) = sign_mat(2,2)*dY; 
new_mat(3,3) = sign_mat(3,3)*dZ; 
new_mat(4,4) = 1; 
new_mat(1,4) = sign_mat(1,4)*newOrigin(1)*dX;
new_mat(2,4) = sign_mat(2,4)*newOrigin(2)*dY;
new_mat(3,4) = sign_mat(3,4)*newOrigin(3)*dZ;

V_new.mat = new_mat;

V_new = spm_write_vol(V_new,Y_new);

close(hw)

end

