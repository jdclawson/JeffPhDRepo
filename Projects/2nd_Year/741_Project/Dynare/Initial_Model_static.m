function [residual, g1, g2, g3] = Initial_Model_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 8, 1);

%
% Model equations
%

T11 = y(4)^params(5);
T45 = y(3)^(-params(3));
lhs =y(1);
rhs =T11;
residual(1)= lhs-rhs;
lhs =y(6);
rhs =params(5)*y(4)^(params(5)-1)-params(2);
residual(2)= lhs-rhs;
lhs =y(7);
rhs =T11*(1-params(5));
residual(3)= lhs-rhs;
lhs =y(2);
rhs =y(7)-y(5)-y(4);
residual(4)= lhs-rhs;
lhs =y(3);
rhs =y(5)*(1+exp(y(8))*params(1))+y(4)*(1+y(6));
residual(5)= lhs-rhs;
lhs =y(2)^(-params(3));
rhs =(1+y(6))*params(4)*T45;
residual(6)= lhs-rhs;
lhs =T11;
rhs =y(5)+y(4)+y(2)+y(3);
residual(7)= lhs-rhs;
lhs =y(8);
rhs =params(1)*(1-params(6))+y(8)*params(6)-x(1);
residual(8)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(8, 8);

  %
  % Jacobian matrix
  %

T64 = getPowerDeriv(y(4),params(5),1);
  g1(1,1)=1;
  g1(1,4)=(-T64);
  g1(2,4)=(-(params(5)*getPowerDeriv(y(4),params(5)-1,1)));
  g1(2,6)=1;
  g1(3,4)=(-((1-params(5))*T64));
  g1(3,7)=1;
  g1(4,2)=1;
  g1(4,4)=1;
  g1(4,5)=1;
  g1(4,7)=(-1);
  g1(5,3)=1;
  g1(5,4)=(-(1+y(6)));
  g1(5,5)=(-(1+exp(y(8))*params(1)));
  g1(5,6)=(-y(4));
  g1(5,8)=(-(y(5)*exp(y(8))*params(1)));
  g1(6,2)=getPowerDeriv(y(2),(-params(3)),1);
  g1(6,3)=(-((1+y(6))*params(4)*getPowerDeriv(y(3),(-params(3)),1)));
  g1(6,6)=(-(params(4)*T45));
  g1(7,2)=(-1);
  g1(7,3)=(-1);
  g1(7,4)=T64-1;
  g1(7,5)=(-1);
  g1(8,8)=1-params(6);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],8,64);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],8,512);
end
end
end
end
