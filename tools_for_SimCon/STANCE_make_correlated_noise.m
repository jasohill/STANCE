function Z = STANCE_make_correlated_noise(nT,CM,PSD_in,dt_in,dt_out)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%
% Generates correlated normalizes noise Z with correlation matrix CM and/or PSD (optional).
%
% Author: Jason E. Hill, Ph.D.
% Date: 25 Jan 2017
%
% nT: number of time points
% dt the sample time step both in and out (optional)
% NOTE: PSD = (abs((fft(X))).^2)/nT;

if nargin < 5
    dt_in = [];
end

if nargin < 6
    dt_out = dt_in;
end

if nargin < 3
    PSD_in = [];
end

if dt_in > dt_out
    error('Output time step must be >/= or the input value.')
end

nR = size(CM,1);
if nR ~= size(CM,2)
   error('Input matrix must be square!');
end

% compute cholesky factor
[cF, p] = chol(CM);
if p > 0
   error('Correlation matrix is not positive semi-definite!')
end

% generate uncorrelated normalized noise
X = STANCE_make_uncorrelated_noise(nT,nR);

% impose power spectral density if it is specified
if ~isempty(PSD_in)
   if (nR ~= size(PSD_in,2)) && (1 ~= size(PSD_in,2)) 
       error('The set of input PSD(s) should be either 1 or the number of ROIs.');
    end  
   % resample PSD to match the number of time points if needed
   if nT ~= size(PSD_in,1)
      for r = 1:size(PSD_in,2) 
         PSD(:,r) = resample(PSD_in(:,r),nT,size(PSD_in,1));
      end
   else
      PSD = PSD_in;
   end

   if ~isempty(dt_in)
      if dt_out ~= dt_in
          % interpolate PSD to ouput sample rate
          for r = 1:size(PSD,2) 
             Vq = ((0:floor(nT/2))*dt_in*length(PSD(:,r)))/(dt_out*(nT))+1;
             PSD_out(:,r) = interp1(PSD(1:ceil(length(PSD(:,r))/2)+1),Vq);
             PSD_out(:,r) = [PSD_out(:,r) fliplr(PSD_out(2:((nT + 1 + mod(nT+1,2))/2-1)),r)]';
             PSD_out(1,r) = 0; 
          end
      else
          PSD_out = PSD;          
      end
   else
       PSD_out = PSD; 
   end

   x = fft(X);
   if size(PSD,2) == 1
       PSD_x = repmat(PSD_out,[1 size(x,2)]);
   else
       PSD_x = PSD_out;
   end
   y = x.*sqrt(PSD_x);
   Y = real(ifft(y));

   % normalize time series to unit variance & zero mean
   Y = (Y-repmat(mean(Y),nT,1))./repmat(std(Y),nT,1);
else
   Y = X;
end

% induce correlation
Z = Y*cF;

% normalize time series to unit variance & zero mean
Z = (Z-repmat(mean(Z),nT,1))./repmat(std(Z),nT,1);

end

