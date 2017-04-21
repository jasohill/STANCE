function bold = balloon(ts_Data, dt, T_scan,  varargin)
% Simulates the balloon model of the hemodynamic response of
% capillaries in the brain.
%   translators:        Xiangyu Liu & Dr. Jason E. Hill (post-doc fellow with CNT at Texas Tech University)
%   date:               5 NOV 2015
%   update:             5 APR 2017
%   source:             neuRosim's balloon() R function
%   source              version: 0.2-12
%   source date:        2015-01-09
%   source authors:     Marijke Welvaert with contributions from Joke Durnez, Beatrijs
%                         Moerkerke, Yves Rosseel and Geert Verdoolaege
%   source maintainer:  Marijke Welvaert <marijke.welvaert@gmail.com>
%
%   for variable options: 'verbose' followed by true or false
%
%   Balloon model parameters (included in param {cell}): 
%                         param.kappa   (default = 2)
%	                      param.tau1    (default = 3)
%	                      param.tauf    (default = 4)
%	                      param.taum    (default = 4)
%	                      param.f1      (default = 1.5)
%	                      param.deltat  (default = 1)
%	                      param.n       (default = 3)
%	                      param.E0      (default = 0.4)
%	                      param.V0      (default = 0.03)
%	                      param.a1      (default = 3.4)
%	                      param.a2      (default = 1.0)
%                         param.tauMTT  (default = 3)
%	                      param.tau     (default = 20)
%	                      param.alpha   (default = 0.4)                           

%% handle variable argument options
nVarargs = length(varargin);
paramOptFlag = 0;

% set default values
param = struct();
param.kappa     = 2;
param.tau1      = 3;
param.tauf      = 4;
param.taum      = 4;
param.f1        = 1.5;
param.deltat    = 1;
param.n         = 3;
param.E0        = 0.4;
param.V0        = 0.03;
param.a1        = 3.4;
param.a2        = 1.0;
param.tauMTT    = 3;
param.tau       = 20;
param.alpha     = 0.4;

verbose = 'false';
if nVarargs > 0
    if mod(nVarargs,2) >0 
        error('optional arguments must come in pairs');
    end
    for i = 1 : 2 : nVarargs
        switch varargin{i}      
            case 'verbose'
                verbose = logical(varargin{i+1});
            case 'param'
                paramOptFlag = 1;
                param = varargin{i+1};
            otherwise
                error('unsupported or unrecognized option');
        end
    end
end

if paramOptFlag == 0
    if verbose
        warning('Default parameter values are used.');
    end
end

%% initializations

param.deltatf = 0;
param.deltatm = param.deltat - param.deltatf;
param.m1 = (param.f1-1)/param.n + 1;

t = dt:dt:T_scan;
it = 1:length(t);
y0 = 0*it;
tspan = [dt T_scan];

param.I = [param.tau1, param.kappa];

[~,I] = ode45(@(t,y) my_inhib(ts_Data,it,y,param.I), tspan, y0);

N = ts_Data - I(end,:)';
F = 1 + conv((param.f1 - 1) * singleGammaHRF(t-param.deltatf, 'FWHM', param.tauf, 'verbose', false), N);
M = 1 + conv((param.m1 - 1) * singleGammaHRF(t-param.deltatm, 'FWHM', param.taum, 'verbose', false), N);
E = param.E0*M./F;
param.balloon = [param.E0, param.tauMTT, param.tau, param.alpha];

y0_2 = ones(2*length(ts_Data),1);
[~,Y] = ode45(@(t,y) de_balloon(E,F,it,y,param.balloon), tspan, y0_2);

V = Y(end,1:length(t));
Q = Y(end,length(t)+1:end);

bold = (param.V0.*(param.a2.*(1 - V) - param.a1.*(1 - Q)));

% normalize
bold = bold/max(bold);

end

% auxillary functions 

function yd1 = my_inhib(ts,it,y,p)
   yd1 = (1/p(1)).*(p(2).*ts(it) - (p(2) + 1).*y(1));
end

function de_ball = de_balloon(E,F,it,y,p)
   dv = (F(it) - y(1).^(1/p(4)))./(p(2)+p(3));
   dq = (1/p(2)).*(F(it).*E(it)./p(1) - y(2)./y(1).*(y(1).^(1/p(4)) + p(3)./(p(2)+p(3)).*(F(it) - y(1).^(1/p(4)))));
   de_ball = [dv;dq];
end


