% 1st Look at my 741 project using dynare, just to get an idea.
% Created by Jeff Clawson

%----------------------------------------------------------------
% 0. Housekeeping
%----------------------------------------------------------------

close all;    % Do not use clear all -- Dynare automatically does that and if you include this command, it would create problems. 

write_latex_dynamic_model;
write_latex_static_model;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y $y$ cy $cy$ co $co$ k $k$ x $x$ r $r$ w $w$ z $z$;
varexo e $e$;

parameters g $g$ delta $\delta$ gamma $\gamma$ beta $\beta$ alpha $\alpha$ rho $rho$;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

% Technology
g  = 0.4;
delta  = 0.012;
gamma = 3;
rho=.9;

% Preferences
beta   = 0.8;
alpha  = 0.35;
sigma = 0.007;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model;
  y = k^alpha;
  r = alpha*k^(alpha-1)-delta;
  w = (1-alpha)*k^(alpha);
  cy = w - x - k;
  co = (1+exp(z(+1))*g)*x+(1+r)*k;
  cy^-gamma = beta*(1+r(+1))*co(+1)^-gamma;
  k^alpha = cy + co + k + x;
  z = (1-rho)*g+ rho*z(-1)-e;
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
  k = 1;
  cy = 0.5;
  co = 0.5;
  z = 0;  
end;

shocks;
var e = sigma^2;
end;

steady;

stoch_simul(hp_filter = 1600, order = 1);

%----------------------------------------------------------------
% 5. Some Results
%----------------------------------------------------------------

