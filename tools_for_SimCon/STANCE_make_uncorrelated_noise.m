function X = STANCE_make_uncorrelated_noise(nT,nR)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%
% Generates totally uncorrelated normalized white noise X.
%
% Jason E. Hill, Ph.D.
% STANCE_make_uncorrelated_noise.m      updated     25 Jan 2017

% create uncorrelated observations

X = randn(nT, nR);

% normalize time series to unit variance & zero mean
X = (X-repmat(mean(X),nT,1))./repmat(std(X),nT,1);

CM = corr(X);
CM(1:(nR + 1):(nR * nR)) = 0;

CMerror = sum((CM(:)).^2)/(nR*(nR-1)/2);

while CMerror > 1e-34
   % compute cholesky factor
   CM = -CM;
   CM(1:(nR + 1):(nR * nR)) = 1;
   cF = chol(CM);

   % induce correlation
   X = X*cF;

   % normalize time series to unit variance & zero mean
   X = (X-repmat(mean(X),nT,1))./repmat(std(X),nT,1);

   CM = corr(X);
   CM(1:(nR + 1):(nR * nR)) = 0;

   CMerror = sum((0 - CM(:)).^2)/(nR*(nR-1)/2);

end

CM(1:(nR + 1):(nR * nR)) = 1;

end