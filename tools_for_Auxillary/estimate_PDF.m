function [pdf,bins] = estimate_PDF(data,bin_size,bin_range)
% estimates the probability density function (PDF) of data over bin_size 
% (default 1.0) over the first data dimension
%
% Author: Jason E. Hill
% Date 27 JAN 2017

if nargin < 2
   bin_size = 1;
end

N = size(data,1);
bindata = round(data/bin_size);

if nargin < 3
   minbin = min(bindata);
   maxbin = max(bindata);
else
   minbin = bin_range(1);
   maxbin = bin_range(2);
end
nbins = maxbin - minbin + 1;

pdf  = zeros(nbins,size(bindata,2));
bins = zeros(nbins,1);
for i = 1:nbins
    member = minbin + i - 1;
    pdf(i,:) = sum(bindata == member)/N;
    bins(i) = member*bin_size;
end

end

