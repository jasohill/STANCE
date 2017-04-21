function delta_x = gamFWHM(k)
%gamFWHM returns the FWHM of the gamma distribution where k = alpha is the
% exponent of x, and theta or beta is taken as one (FWHM ~ theta)
% NOTE: must have k > 1 for a positive and meaningful result.
%
% Jason E. Hill, Post-Doc with the CNT at Texas Tech.
% gamFWHM   date     25 MAR 2017

delta_x = ((k-1)*(lambertw(-2^(1/(1-k))/exp(1))-lambertw(-1,-2^(1/(1-k))/exp(1))));

end

