% Basic RBC Model as described in Chapter 1 of Cooley, 
% Frontiers of Business Cycle Research
% 
% Details:
%   1) Computation in levels, 1st order approximation.
%   2) Calculates also by perturbation investment (i) and labor 
%   productivity (y_l).
%   3) Calibration from Cooley and Prescott.
%
% Jesus Fernandez-Villaverde
% Philadelphia, March 3, 2005

% Edited by Boragan Aruoba

%----------------------------------------------------------------
% 0. Housekeeping
%----------------------------------------------------------------

close all;    % Do not use clear all -- Dynare automatically does that and if you include this command, it would create problems. 

write_latex_dynamic_model;
write_latex_static_model;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y $y$ c $c$ k $k$ x $x$ h $h$ p $p$ z $z$;
varexo e $e$;

parameters theta $\theta$ delta $\delta$ rho $\rho$ sigma $\sigma$ gamma $\gamma$ beta $\beta$ alpha $\alpha$ eta $\eta$;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

% Technology
theta  = 0.4;
delta  = 0.012;
rho    = 0.95;
sigma  = 0.007; 
gamma  = (1.0156^0.25)-1;

% Preferences
beta   = 0.987;
alpha  = 0.64;
eta    = (1.012^0.25)-1;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model;
  ((1+gamma)/c) = beta*(1/c(+1))*(1+theta*exp(z(+1))*(k^(theta-1))*(h(+1)^(1-theta))-delta);
  alpha*c/((1-alpha)*(1-h)) = (1-theta)*exp(z)*(k(-1)^theta)*(h^(-theta));
  c+x = y;
  y = exp(z)*(k(-1)^theta)*(h^(1-theta));
  x = (1+gamma)*(1+eta)*k-(1-delta)*k(-1);
  p = y/h;
  z = rho*z(-1)+e;
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
  k = 24;
  c = 1.33;
  h = 0.31;
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

statistic1 = 100*sqrt(diag(oo_.var(1:6,1:6)))./oo_.mean(1:6);
dyntable('Relative standard deviations in %',strvcat('VARIABLE','REL. S.D.'),M_.endo_names(1:6,:),statistic1,10,8,4);
