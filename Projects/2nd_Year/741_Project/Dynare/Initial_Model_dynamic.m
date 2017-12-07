function [residual, g1, g2, g3] = Initial_Model_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(8, 1);
T11 = y(5)^params(5);
T48 = y(10)^(-params(3));
lhs =y(2);
rhs =T11;
residual(1)= lhs-rhs;
lhs =y(7);
rhs =params(5)*y(5)^(params(5)-1)-params(2);
residual(2)= lhs-rhs;
lhs =y(8);
rhs =T11*(1-params(5));
residual(3)= lhs-rhs;
lhs =y(3);
rhs =y(8)-y(6)-y(5);
residual(4)= lhs-rhs;
lhs =y(4);
rhs =y(6)*(1+exp(y(12))*params(1))+y(5)*(1+y(7));
residual(5)= lhs-rhs;
lhs =y(3)^(-params(3));
rhs =params(4)*(1+y(11))*T48;
residual(6)= lhs-rhs;
lhs =T11;
rhs =y(6)+y(5)+y(3)+y(4);
residual(7)= lhs-rhs;
lhs =y(9);
rhs =params(1)*(1-params(6))+params(6)*y(1)-x(it_, 1);
residual(8)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(8, 13);

  %
  % Jacobian matrix
  %

T69 = getPowerDeriv(y(5),params(5),1);
  g1(1,2)=1;
  g1(1,5)=(-T69);
  g1(2,5)=(-(params(5)*getPowerDeriv(y(5),params(5)-1,1)));
  g1(2,7)=1;
  g1(3,5)=(-((1-params(5))*T69));
  g1(3,8)=1;
  g1(4,3)=1;
  g1(4,5)=1;
  g1(4,6)=1;
  g1(4,8)=(-1);
  g1(5,4)=1;
  g1(5,5)=(-(1+y(7)));
  g1(5,6)=(-(1+exp(y(12))*params(1)));
  g1(5,7)=(-y(5));
  g1(5,12)=(-(y(6)*exp(y(12))*params(1)));
  g1(6,3)=getPowerDeriv(y(3),(-params(3)),1);
  g1(6,10)=(-(params(4)*(1+y(11))*getPowerDeriv(y(10),(-params(3)),1)));
  g1(6,11)=(-(params(4)*T48));
  g1(7,3)=(-1);
  g1(7,4)=(-1);
  g1(7,5)=T69-1;
  g1(7,6)=(-1);
  g1(8,1)=(-params(6));
  g1(8,9)=1;
  g1(8,13)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],8,169);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],8,2197);
end
end
end
end
