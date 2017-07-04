function [V_signal_out,Y_signal_out] = STANCE_EPI_signal(fn_tissues,T2star,scan,noiseMap,physio4D,motion,approximation,gzipFlag,seed)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% A signal equation solver for EPI sequences.
%
%    fn_tissue:     filename of tissue priors for target space
%    T2star:        T2star data of source (3D or 4D)
%    noiseMap:      a 3D mask allowing for spatial variation in noise via
%                   tissue priors
%    motion:        [translation_x,translation_y,translation_z,rotation_x,...
%                      rotation_y,rotation_z] where rotation_z is axial in-plane rotation 
%    approximation: flag for computation approximation (less accurate but faster) 
%    gzipFlag:      flag for saving in the archive file format *.gz
%
% author: Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
% STANCE_EPI_signal.m      updated     23 MAR 2017
%
%------------------------------------------------------------------------
% 20 Anatomical Models of 20 Normal Brains from
% "An new improved version of the realistic digital brain phantom" 
% - Berengere Aubert-Broche, Alan C. Evans, & Louis Collins
%   NeuroImage 32 (2006) 138 - 145
%------------------------------------------------------------------------
% NAME              relative PD     T1[ms]      T2[ms]      T2*[ms]
% Background*       0    5.0e-4          7.3e-4      0.3         0.001
% CSF               1    1.0             2569.0      329.0       58.0
% Grey Matter       2    0.86            833.0       83.0        69.0
% White Matter      3    0.77            500.0       70.0        61.0
% Fat               4    1.0             350.0       70.0        58.0
% Muscle            5    1.0             900.0       47.0       *40.0
% Muscle/Skin       6    1.0            *1800.0     *140.0      *48.0
% Skull*            7    0.28            2580.0      314.8       50.0
% Vessels*          8    0.87            900.0       180.0       58.0
% Connective (fat2) 9    0.77            500.0       70.0        61.0
% Dura Mater       10    1.0             2569.0      329.0       58.0
% Bone Marrow      11    1.0            *2500.0      70.0        61.0

% handle variable input arguments
if nargin<9
    % do nothing, carry on
else
    if ~isempty(seed)
        rng(seed,'twister');
    end
end
if nargin < 8
    gzipFlag = false;
else
    if isempty(gzipFlag)
        gzipFlag = false;
    end    
end
if nargin < 7
    approximation = false;
else
    if isempty(approximation)
        approximation = false;
    end        
end
if nargin < 6
    motion = [];
end
if nargin < 5
    physio4D = [];
end
if nargin < 4
    noiseMap = [];
end
if nargin < 3
    scan = [];
end
if nargin < 3 || isempty(scan)
    scan.voxel.size    = [3 3 3]; % [mm] voxel size dimensions
    scan.voxel.matrix  = [64 64 NaN]; % size of matrix (number of Z slices may be autodetected)
    scan.voxel.spacing = [0 0 0.2*scan.voxel.size(3)]; % assume 20% Z spacing
    scan.tiltAngle     = 15;   % [degrees] tilt angle from AC-PC line
    scan.TR            = 2000; % [ms] repetion time
    scan.TE            = 30;   % [ms] echo time
    scan.ES            = 0.51; % [ms] echo spacing
    scan.FA            = 90;   % [degrees] flip angle 
    scan.BW            = 2232; % [Hz/Px] bandwidth
    scan.order         = 'SD'; % SD = sequential descending order
    scan.B0             = 3.0; % [T] the main magnet field strength        
    scan.KM0           = 2225; % fit to data with max of 909 at 3T and FA = 90 degrees
    scan.noise_method  = 'sigma';
    scan.noise         = 0;   
    scan.attenuation   = 0;
    scan.acceleration  = 1;
end
if ~isfield(scan,'voxel')
    scan.voxel.size    = [3 3 3]; % [mm] voxel size dimensions
    scan.voxel.matrix  = [64 64 NaN];  % size of matrix (number of Z slices may be autodetected)
    scan.voxel.spacing = [0 0 0.2*scan.voxel.size(3)]; % assume 20% Z spacing
end
if ~isfield(scan.voxel,'size')
    scan.voxel.size    = [3 3 3]; % [mm] voxel size dimensions
end
if ~isfield(scan.voxel,'matrix')
    scan.voxel.matrix  = [64 64 NaN]; % size of matrix (number of Z slices may be autodetected)
end
if ~isfield(scan.voxel,'spacing')
    scan.voxel.spacing = [0 0 0.2*scan.voxel.size(3)]; % assume 20% Z spacing
end
if ~isfield(scan,'tiltAngle')
    scan.tiltAngle     = 15; % [degrees] 
end
if ~isfield(scan,'TR')
    scan.TR            = 2000; % [ms]
end
if ~isfield(scan,'TE')
    scan.TE            = 30;   % [ms]
end
if ~isfield(scan,'ES')
    scan.ES            = 0.51; % [ms] echo spacing
end
if ~isfield(scan,'FA')
    scan.FA            = 90;   % degrees 
end
if ~isfield(scan,'BW') 
    scan.BW            = 2232; % [Hz/Px] bandwidth
end
if ~isfield(scan,'order') 
    scan.order         = 'SD'; % SD = sequential descending order
end
if ~isfield(scan,'B0') 
    scan.B0             = 3.0; % [T] the main magnet field strength    
end
if ~isfield(scan,'KM0') 
    scan.KM0           = 2225; % fit to data with max of 909 at 3T and FA = 90 degrees
end
if ~isfield(scan,'noise_method') 
    scan.noise_method  = 'sigma';
end
if ~isfield(scan,'noise_method') 
    scan.noise_method  = 'sigma';
end
if ~isfield(scan,'noise') 
    scan.noise         = 0;   
end
if ~isfield(scan,'attenuation')   
    scan.attenuation   = 0;
end
if ~isfield(scan,'acceleration')   
    scan.acceleration   = 1;
end

TR          = scan.TR;
TE          = scan.TE;
FA          = scan.FA;
KM0         = scan.KM0;
method      = scan.noise_method;
noise       = scan.noise;
attenuation = scan.attenuation;
sliceOrder  = scan.order;
acceleration  = scan.acceleration;

S = sind(FA);
C = cosd(FA);

% generate T1 relaxation time and proton density (PD) baseline maps
[V_tissues,Y_tissues] = STANCE_load_volume(fn_tissues);
[~,T1_Map] = STANCE_make_parameter_map(fn_tissues,'T1');
[~,PD_Map] = STANCE_make_parameter_map(fn_tissues,'PD');

if ~exist('STANCE.mat') %#ok<*EXIST>
    if ~exist('../STANCE.mat')
       load('../../STANCE.mat');
    else
       load('../STANCE.mat');
    end
else
    load('STANCE.mat');
end


sz = size(T2star);
if length(sz) < 4
    sz(4) = 1;
end
sliceTiming = make_slice_timing(sliceOrder,sz(3),acceleration);
N_ADCs = length(sliceTiming); 
if isempty(motion)
    % do nothing
elseif size(motion,1) == 1 && sz(4) == 1
        
   if size(motion,2) > 3
      rotation_x = motion(4);
      rotation_y = motion(5);
      rotation_z = motion(6); 
   end
   translation_x = motion(1);
   translation_y = motion(2);
   translation_z = motion(3);
       
elseif sz(4) > 1 
    sz_motion = size(motion);
    motion(sz_motion(1)+1,:) = motion(sz_motion(1),:);
    if sz_motion(2) > 3
       rotation_x = 0.5*(motion(1:sz_motion(1),4)+motion(2:sz_motion(1)+1,4));
       rotation_y = 0.5*(motion(1:sz_motion(1),5)+motion(2:sz_motion(1)+1,5));
       rotation_z = motion(:,6); 
       dtheta_x = (pi/360)*(motion(2:sz_motion(1)+1,4) - motion(1:sz_motion(1),4));
       dtheta_y = (pi/360)*(motion(2:sz_motion(1)+1,5) - motion(1:sz_motion(1),5));
    end
    translation_x = 0.5*(motion(1:sz_motion(1),1)+motion(2:sz_motion(1)+1,1));
    translation_y = 0.5*(motion(1:sz_motion(1),2)+motion(2:sz_motion(1)+1,2));
    translation_z = 0.5*(motion(1:sz_motion(1),3)+motion(2:sz_motion(1)+1,3));
    
else
   warning('Mismatch between motion and data dimensions!')
end
       

Y_signal_out = 0*T2star;
hw = waitbar(0,'Computing EPI signal timeseries ...');
for t = 1:sz(4)
    if sz(4) == 1
        sz2 = size(T2star);
        if length(sz2) < 4
            T2star_Map = T2star;
        else        
            T2star_Map = squeeze(T2star(:,:,:,1));    
        end
    else
        T2star_Map = squeeze(T2star(:,:,:,t));  
    end

    if ~approximation
        E1 = exp(-TR./T1_Map);    
        Y_signal = KM0.*PD_Map.*(S.*(1-E1)./(1-C.*E1)).*exp(-TE/T2star_Map);
    else
        a1 = TR./833.0; % for gray matter average
        a2 = TE./69.0;  % for gray matter average
    
        R1 = TR./T1_Map;
        R2 = TE./T2star_Map;    

        Ea1 = exp(-a1);
        Ea2 = exp(-a2);
    
        D = 1.0/(1 - C*Ea1);
        A1 = (1-Ea1)*D;
        B1 = Ea1*(1-C)*D^2;
        C1 = -0.5*B1*(1+C*Ea1)*D;
    
        Y_signal = KM0.*PD_Map.*S.*(A1 + B1*(R1-a1) + C1*(R1-a1).^2).*Ea2.*(1-(R2-a2)+0.5*(R2-a2).^2);

    end

    Y_signal(isnan(Y_signal(:))) = 0;

    if isempty(physio4D)
        %do nothing more
    else 
        Y_signal = Y_signal.*(1 + physio4D(:,:,:,t));
    end

    if t == 1
        V_signal = V_tissues(1);
    end


%% add motion
    if isempty(motion)
        % do nothing
    elseif size(motion,1) == 1 && sz(4) == 1
        
       if size(motion,2) > 3  
%           tic
           if ~ismatrix(Y_signal)
               for x = 1:size(Y_signal,1)
                  I_x_slice = squeeze(Y_signal(x,:,:));
                  Y_signal(x,:,:) = imrotate(I_x_slice,rotation_x,'bilinear','crop');
               end               
               for y = 1:size(Y_signal,2)
                  I_y_slice = squeeze(Y_signal(:,y,:));
                  Y_signal(:,y,:) = imrotate(I_y_slice,rotation_y,'bilinear','crop');
               end
               for z = 1:size(Y_signal,3)
                  I_z_slice = Y_signal(:,:,z,t);
                  Y_signal(:,:,z) = imrotate(I_z_slice,rotation_z,'bilinear','crop');
               end
           else
               Y_signal = imrotate(Y_signal,rotation_z,'bilinear','crop');           
           end
%           toc
% Note: MATLAB 2017a+ users can use the following:           
%            tic
%            if ndims(Y_signal)>2
%                for x = 1:size(Y_signal,1)
%                   Y_signal = imrotate3(Y_signal,rotation_x,[1 0 0],'linear','crop','FillValues',0);
%                end
%                for y = 1:size(Y_signal,2)
%                   Y_signal = imrotate3(Y_signal,rotation_y,[0 1 0],'linear','crop','FillValues',0);
%                end
%                for z = 1:size(Y_signal,3)
%                   Y_signal = imrotate3(Y_signal,rotation_z,[0 0 1],'linear','crop','FillValues',0);
%                end
%            else
%                Y_signal = imrotate(Y_signal,rotation_z,'bilinear','crop');           
%            end
%            toc           
       end
       Y_signal = imtranslate(Y_signal,[translation_x,translation_y,translation_z]);
       
    elseif sz(4) > 1
    
       if size(motion,2) > 3
           rotation_x = motion(:,4);
           rotation_y = motion(:,5);
           rotation_z = motion(:,6);  
%           tic
               for x = 1:size(Y_signal,1)
                  I_x_slice = squeeze(Y_signal(x,:,:));
                  % reduce magnitude due to average spin-history over TR:
                  %    (1-(b/t)*\theta/N_ADCs) where b = base, t = slice thickness
                  %    and using tan(\theta) ~ \theta for small angles
                  I_x_slice = I_x_slice*(1-(0.5/N_ADCs)*dtheta_x(t)*scan.voxel.size(1)*scan.voxel.matrix(1)/scan.voxel.size(3));
                  % add rigid body rotation
                  Y_signal(x,:,:) = imrotate(I_x_slice,rotation_x(t),'bilinear','crop');
               end               
               for y = 1:size(Y_signal,2)
                  I_y_slice = squeeze(Y_signal(:,y,:));
                  % reduce magnitude due to average spin-history over TR:
                  %    (1-(b/s)*\theta/N_ADCs) where b = base, s = slice thickness
                  %    and using tan(\theta) ~ \theta for small angles
                  I_y_slice = I_y_slice*(1-(0.5/N_ADCs)*dtheta_y(t)*scan.voxel.size(2)*scan.voxel.matrix(2)/scan.voxel.size(3));
                  % add rigid body rotation                  
                  Y_signal(:,y,:) = imrotate(I_y_slice,rotation_y(t),'bilinear','crop');
               end
               for z = 1:size(Y_signal,3)
                  [~,col] = find(sliceTiming == z);
                  I_z_slice = squeeze(Y_signal(:,:,z));
                  % add spin-history magnitude reduction slice-by-slice
                  dz0 = motion(t+1,3)*((col-1)/N_ADCs)-motion(t,3)*(1-(col-1)/N_ADCs);
                  dzTA = motion(t+1,3)*((col)/N_ADCs)-motion(t,3)*(1-(col)/N_ADCs);
                  dz = 0.5*(dzTA-dz0);
                  I_z_slice = I_z_slice*(1-dz/scan.voxel.size(3));
                  % add rigid body rotation 
                  Y_signal(:,:,z) = imrotate(I_z_slice,rotation_z(t)*(1-(col-1)/N_ADCs)+rotation_z(t+1)*((col-1)/N_ADCs),'bilinear','crop');
               end
%           toc
% Note: MATLAB 2017a+ users can use the following:           
%            tic
%            for t = 1:size(motion,1)
%                for x = 1:size(Y_signal,1)
%                   Y_signal(:,:,:,t) = imrotate3(squeeze(Y_signal(:,:,:,t)),rotation_x,[1 0 0],'linear','crop','FillValues',0);
%                end
%                for y = 1:size(Y_signal,2)
%                   Y_signal(:,:,:,t) = imrotate3(squeeze(Y_signal(:,:,:,t)),rotation_y,[0 1 0],'linear','crop','FillValues',0);
%                end
%                for z = 1:size(Y_signal,3)
%                   Y_signal(:,:,:,t) = imrotate3(squeeze(Y_signal(:,:,:,t)),rotation_z,[0 0 1],'linear','crop','FillValues',0);
%                end
%            end
%            toc           
       end
       Y_signal = imtranslate(Y_signal,[translation_x(t),translation_y(t),translation_z(t)]);        
    else
       warning('Mismatch between motion and data dimensions!')
    end


%% add noise

    if noise > 0
        if strcmp(method, 'percent');
            S_max = max(Y_signal(:));
            sigma = S_max.*noise/100;
        else
            sigma = noise;
        end
        [M, N, L] = size(Y_signal);
        if isempty(noiseMap)
            noiseVolume = sigma*randn(M,N,L) + 1j*sigma*randn(M,N,L);
        else
            noiseVolume = (sigma*randn(M,N,L) + 1j*sigma*randn(M,N,L)).*noiseMap;
        end
        Y_signal = abs(Y_signal + noiseVolume);
    end

%% add attentuation

    if attenuation > 0
        stopFlag = false;
        brainMask0 = (1 - Y_tissues(:,:,:,1));
        brainMask = brainMask0;
        lastSum = sum(brainMask(:));
        for n = 1:40
            w = 2*n+1;
            if ~stopFlag
                update = imerode(brainMask0,ones(w,w,w));
                if sum(update(:)) == lastSum
                    stopFlag = true;
                else
                   brainMask = brainMask + update;
                   lastSum = sum(update(:));
                end
            end
        end
        attenuationMask = exp(-brainMask./attenuation);
        Y_signal = attenuationMask.*Y_signal;
    end    

    
    % save volume for time t
    if t == 1
        V_signal_out = V_signal;
    end
    Y_signal_out(:,:,:,t) = Y_signal;
    waitbar(t/sz(4),hw)    
    
end


%% save results
fullpathstr = STANCE_genpath;
if ~exist(fullpathstr,'dir')
   STANCE_new_session;
end

%    V_signal.fname = sprintf('%s_%s_%s_%.3d%s','EPI_BOLD_',fullpathstr(end-14:end-11),fullpathstr(end-2:end),t,'.nii');
%    fprintf('Writing %s\n',V_signal.fname);
    
V_signal.fname = [STANCE_genpath,'\EPI_BOLD_',fullpathstr(end-14:end-11),'_',fullpathstr(end-2:end),'.nii']; 
fprintf('o Writing %s\n',V_signal.fname);
V_signal.private.dat.fname = V_signal.fname;
if length(V_signal.private.dat.dim) > 3
    V_signal.private.dat.dim = V_signal.private.dat.dim(1:3);
end

cd(fullpathstr);

V_signal = spm_write_vol(V_signal,single(Y_signal));

nii = load_nii(V_signal.fname);

nii.img = Y_signal_out;
   
% Update header    
nii.hdr.dime.dim(1) = 4;
nii.hdr.dime.dim(5) = sz(4);
nii.hdr.dime.pixdim(5) = scan.TR;
max_out = max(Y_signal_out(:));
X = ['The maximum intensity of the simulated signal: ',num2str(max_out)];
disp(X)

nii.hdr.dime.glmax = max(Y_signal_out(:));
nii.hdr.dime.glmin = 0.0;

sliceOrder = scan.order;
nii.hdr.hist.slice_code = 0;  % default value
    switch sliceOrder
        case 'SA'
            nii.hdr.hist.slice_code = 1;        
        case 'SD'
            nii.hdr.hist.slice_code = 2; 
        case 'IA'
            nii.hdr.hist.slice_code = 3;           
        case 'ID'
            nii.hdr.hist.slice_code = 4; 
        case 'IA2'
            nii.hdr.hist.slice_code = 5;            
        case 'ID2'
            nii.hdr.hist.slice_code = 6;            
    end

    
delete(V_signal.fname);
save_nii(nii,V_signal.fname);
   
if gzipFlag
    gzip(V_signal.fname);
    delete(V_signal.fname);
end
cd(STANCEroot);

Y_signal_out = squeeze(Y_signal_out);

close(hw)

end

