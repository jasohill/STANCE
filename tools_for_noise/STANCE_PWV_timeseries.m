function v_ts = STANCE_PWV_timeseries(tC_ts,PW_width,PW_width_sigma,v_0,v_constant_fraction,seed)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
%
% Generates a pulse wave velocity time-series characterized by 
% cardiac impulse time-series tR_ts with pulse wave width PW_width
% width average blood flow v_0 and constant fraction v_constant_fraction. 
% Output is normalized to average amplitude.
%
%   author:         Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
%   updated:        24 MAR 2017

if nargin<6
    % do nothing, carry on
else
    if ~isempty(seed)
        rng(seed,'twister');
    end
end
if nargin < 5
    v_constant_fraction = 0.1;
end
if isempty(v_constant_fraction)
    v_constant_fraction = 0.1;
end
if nargin < 4
    v_0 = 3.0; % mm/s
end
if isempty(v_0)
    v_0 = 3.0; % mm/s
end
if nargin < 3
    PW_width_sigma = 0.02;
end
if isempty(PW_width_sigma)
    PW_width_sigma = 0.02;
end

if nargin < 2
    PW_width = 0.5;
end
if isempty(PW_width)
    PW_width = 0.5;
end

% extract the time intervals from the cardiac impulse timeseries 
Nt = length(tC_ts.Data);
dt = tC_ts.Time(2)-tC_ts.Time(1);
CT_indeces = find(tC_ts.Data == 1);
Nbeats = length(CT_indeces);
CTI = (Nt*dt)/Nbeats;
CIs = CT_indeces(2:Nbeats) - CT_indeces(1:Nbeats-1);
CIs(Nbeats) = CIs(Nbeats-1);
CIs(2:Nbeats-1) = round(0.5*(CIs(1:Nbeats-2) + CIs(2:Nbeats-1)));
CIs = 2.*round(0.5.*CIs)+1; %ensure odd durations 

% A single pulse wave is represented with a functional form:
% v(phi) = |(cos(2*pi*phi)|
% As detailed in the description of the aorta pulse wave velocity in
% Section 5 of "Numerical modelling of pulse wave propagation in the cardiovascular system: development, validation and clinical applications" 
% the Ph.D. Thesis of J.A. Arimon, The University of London (2006).
% As modelled empirically from data shown in the report
% "2. Cerebral blood flow and Metabolism" in 02-Cerebral_blood_flow.pdf
z_ts_Data = v_constant_fraction*v_0*ones(1,Nt);
As = 0.01.*randn(Nbeats,1) + 2*(1-v_constant_fraction)*pi*v_0*(CIs*dt)./CTI;
PWs = PW_width + PW_width_sigma.*randn(Nbeats,1);
if (CT_indeces(1) - 0.5*(round(PWs(1)*CIs(1))+1)) < 1
    start_i = 1 + 0.5.*(round(PWs(1)*CIs(1))+1)-CT_indeces(1);
    phi = start_i:round(PWs(1)*CIs(1));
    z_phi = As(1).*cos((pi/PWs(1)).*(phi - 0.5.*(round(PWs(1)*CIs(1))+1))./CIs(1));
    z_ts_Data(1:length(phi)) = z_phi + z_ts_Data(1:length(phi)); 
else
    phi = 1:round(PWs(1)*CIs(1));
    z_phi = As(1).*cos((pi/PWs(1)).*(phi - 0.5.*(round(PWs(1)*CIs(1))+1))./CIs(1));
    z_ts_Data(round(CT_indeces(1) - 0.5*(round(PWs(1)*CIs(1))-1)):round(CT_indeces(1) + 0.5*(round(PWs(1)*CIs(1))-1))) = z_phi + z_ts_Data(round(CT_indeces(1) - 0.5*(round(PWs(1)*CIs(1))-1)):round(CT_indeces(1) + 0.5*(round(PWs(1)*CIs(1))-1))); 
end 
for n = 2:(Nbeats-1)
    start_index_ts = round(CT_indeces(n) - 0.5*(round(PWs(n)*CIs(n))-1));
    start_index_phi = 1;
    if start_index_ts < 1
        start_index_phi = start_index_phi + (1 - start_index);
        start_index_ts = round(CT_indeces(n) + 0.5*(round(PWs(n)*CIs(n))-1)) - round(PWs(n)*CIs(n)) + start_index_phi;
    end
    if round(CT_indeces(n) + 0.5*(round(PWs(n)*CIs(n))-1)) > length(z_ts_Data)
        end_index_ts = length(z_ts_Data);
    else
        end_index_ts = round(CT_indeces(n) + 0.5*(round(PWs(n)*CIs(n))-1));
    end
    end_index_phi = start_index_phi + (end_index_ts-start_index_ts);
    phi = start_index_phi:end_index_phi;
    z_phi = As(n).*cos((pi/PWs(n)).*(phi - 0.5.*(round(PWs(n)*CIs(n))+1))./CIs(n));
    z_ts_Data(start_index_ts:end_index_ts) = z_phi + z_ts_Data(start_index_ts:end_index_ts); 
end
if (CT_indeces(Nbeats) + 0.5*(CIs(Nbeats)+1)) > Nt
    end_i = round(Nt-0.5.*PWs(n).*(CIs(Nbeats)+1)-CT_indeces(Nbeats));
    phi = 1:end_i;
    z_phi = As(Nbeats).*cos((pi/PWs(Nbeats)).*(phi - 0.5.*(round(PWs(Nbeats)*CIs(Nbeats))+1))./CIs(Nbeats));
    z_ts_Data(round(CT_indeces(Nbeats) - 0.5*(round(PWs(Nbeats)*CIs(Nbeats))+1)):end_i) = z_phi + z_ts_Data(round(CT_indeces(Nbeats) - 0.5*(round(PWs(Nbeats)*CIs(Nbeats))-1)):end_i); 
else
    phi = 1:round(PWs(Nbeats)*CIs(Nbeats));
    z_phi = As(Nbeats).*cos((pi/PWs(Nbeats)).*(phi - 0.5.*(round(PWs(Nbeats)*CIs(Nbeats))+1))./CIs(Nbeats));
    z_ts_Data((CT_indeces(Nbeats) - round(0.5*(round(PWs(Nbeats)*CIs(Nbeats))-1))):(CT_indeces(Nbeats) + 0.5*(round(PWs(Nbeats)*CIs(Nbeats))-1))) = z_phi + z_ts_Data((CT_indeces(Nbeats) - 0.5*(round(PWs(Nbeats)*CIs(Nbeats))-1)):(CT_indeces(Nbeats) + round(0.5*(round(PWs(Nbeats)*CIs(Nbeats))-1)))); 
end 

ts_time = squeeze(dt:dt:(Nt*dt));
v_ts = timeseries(squeeze(z_ts_Data),ts_time);
v_ts.TimeInfo.Units = 'seconds';
v_ts = setuniformtime(v_ts,'StartTime',dt);
v_ts = setuniformtime(v_ts,'Interval',dt);

