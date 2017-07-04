function z_ts = make_pulse_timeseries(dt,Nt,TI,TI_sigma,p,A_sigma,SighTI,SighTI_sigma,impulse_seed,dz_seed)
% Generates a timeseries of pulsation events characterized by 
% interval time TI between pulses with standard deviation TI_sigma
% normal relaxed position is baseline z = 0. 
% Shape of pulses are half cosines raised to the p-power
% Output is normalized to average amplitude.
% Defaults assume a typical respiratory pulse.
%
%   author:         Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
%   updated:        5 APR 2017

if nargin < 10
    dz_seed = [];
end
if nargin < 9
    impulse_seed = [];
end
if nargin < 8
    SighTI_sigma = 15;
end
if isempty(SighTI_sigma)
    SighTI_sigma = 15;
end    
if nargin < 7
    SighTI = 400;
end
if isempty(SighTI)
    SighTI = 400;
end    
if nargin < 6
    A_sigma = 0.005;
end
if isempty(A_sigma)
    A_sigma = 0.005;
end
if nargin < 5
    p = 4.0;
end
if isempty(p)
    p = 4.0;
end
if nargin < 4
    TI_sigma = 0.25;
end
if isempty(TI_sigma)
    TI_sigma = 0.25;
end
if nargin < 3
    TI = 4.0;
end
if isempty(TI)
    TI = 4.0;
end

% make the impulse timeseries
impulse_ts = make_impulse_timeseries(dt,Nt,TI,TI_sigma,impulse_seed);
indices = find(impulse_ts.Data == 1);
N = length(indices);
% calculate pulse intervals
PIs = indices(2:N) - indices(1:N-1);
PIs(N) = PIs(N-1);
PIs(2:N-1) = round(0.5*(PIs(1:N-2) + PIs(2:N-1)));
PIs = 2.*round(0.5.*PIs)+1; %ensure odd durations 

z_ts_Data = zeros(1,Nt);

if ~isempty(dz_seed)
    rng(dz_seed,'twister')
end
As = A_sigma.*randn(N,1) + 1.0.*(PIs*dt)./TI;
if (indices(1) - 0.5*(PIs(1)+1)) < 1
    start_i = 1 + 0.5*(PIs(1)+1)-indices(1);
    phi = start_i:PIs(1);
    z_phi = As(1).*cos(pi.*(phi - 0.5.*(PIs(1)+1))./PIs(1)).^p;
    z_ts_Data(1:length(phi)) = z_phi + z_ts_Data(1:length(phi)); 
else
    phi = 1:PIs(1);
    z_phi = As(1).*cos(pi.*(phi - 0.5.*(PIs(1)+1))./PIs(1)).^p;
    z_ts_Data((indices(1) - 0.5*(PIs(1)-1)):(indices(1) + 0.5*(PIs(1)-1))) = z_phi + z_ts_Data((indices(1) - 0.5*(PIs(1)-1)):(indices(1) + 0.5*(PIs(1)-1))); 
end 
for n = 2:(N-1)
    phi = 1:PIs(n);
    z_phi = As(n).*cos(pi.*(phi - 0.5.*(PIs(n)+1))./PIs(n)).^p;
    z_ts_Data((indices(n) - 0.5*(PIs(n)-1)):(indices(n) + 0.5*(PIs(n)-1))) = z_phi + z_ts_Data((indices(n) - 0.5*(PIs(n)-1)):(indices(n) + 0.5*(PIs(n)-1))) ; 
end
if (indices(N) + 0.5*(PIs(N)+1)) > Nt
    end_i = Nt-0.5.*(PIs(N)+1)-indices(N);
    phi = 1:end_i;
    z_phi = As(N).*cos(pi.*(phi - 0.5.*(PIs(N)+1))./PIs(N)).^p;
    z_ts_Data((indices(N) - 0.5*(PIs(N)+1)):end_i) = z_phi + z_ts_Data((indices(N) - 0.5*(PIs(N)-1)):end_i); 
else
    phi = 1:PIs(N);
    z_phi = As(N).*cos(pi.*(phi - 0.5.*(PIs(N)+1))./PIs(N)).^p;
    z_ts_Data((indices(N) - 0.5*(PIs(N)-1)):(indices(N) + 0.5*(PIs(N)-1))) = z_phi + z_ts_Data((indices(N) - 0.5*(PIs(N)-1)):(indices(N) + 0.5.*(PIs(N)-1))); 
end 

ts_time = squeeze(dt:dt:(Nt*dt));
z_ts = timeseries(squeeze(z_ts_Data),ts_time);
z_ts.TimeInfo.Units = 'seconds';
z_ts = setuniformtime(z_ts,'StartTime',dt);
z_ts = setuniformtime(z_ts,'Interval',dt);

