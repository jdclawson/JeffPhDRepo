%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'Initial_Model';
M_.dynare_version = '4.5.1';
oo_.dynare_version = '4.5.1';
options_.dynare_version = '4.5.1';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('Initial_Model.log');
M_.exo_names = 'e';
M_.exo_names_tex = 'e';
M_.exo_names_long = 'e';
M_.endo_names = 'y';
M_.endo_names_tex = 'y';
M_.endo_names_long = 'y';
M_.endo_names = char(M_.endo_names, 'cy');
M_.endo_names_tex = char(M_.endo_names_tex, 'cy');
M_.endo_names_long = char(M_.endo_names_long, 'cy');
M_.endo_names = char(M_.endo_names, 'co');
M_.endo_names_tex = char(M_.endo_names_tex, 'co');
M_.endo_names_long = char(M_.endo_names_long, 'co');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'x');
M_.endo_names_tex = char(M_.endo_names_tex, 'x');
M_.endo_names_long = char(M_.endo_names_long, 'x');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names_long = char(M_.endo_names_long, 'z');
M_.endo_partitions = struct();
M_.param_names = 'g';
M_.param_names_tex = 'g';
M_.param_names_long = 'g';
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, '\delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'gamma');
M_.param_names_tex = char(M_.param_names_tex, '\gamma');
M_.param_names_long = char(M_.param_names_long, 'gamma');
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, '\beta');
M_.param_names_long = char(M_.param_names_long, 'beta');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, '\alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names_long = char(M_.param_names_long, 'rho');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 8;
M_.param_nbr = 6;
M_.orig_endo_nbr = 8;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
erase_compiled_function('Initial_Model_static');
erase_compiled_function('Initial_Model_dynamic');
M_.orig_eq_nbr = 8;
M_.eq_nbr = 8;
M_.ramsey_eq_nbr = 0;
M_.lead_lag_incidence = [
 0 2 0;
 0 3 0;
 0 4 10;
 0 5 0;
 0 6 0;
 0 7 11;
 0 8 0;
 1 9 12;]';
M_.nstatic = 5;
M_.nfwrd   = 2;
M_.npred   = 0;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 1;
M_.ndynamic   = 3;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(8, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(6, 1);
M_.NNZDerivatives = [25; -1; -1];
close all;    
M_.params( 1 ) = 0.4;
g = M_.params( 1 );
M_.params( 2 ) = 0.012;
delta = M_.params( 2 );
M_.params( 3 ) = 3;
gamma = M_.params( 3 );
M_.params( 6 ) = .9;
rho = M_.params( 6 );
M_.params( 4 ) = 0.8;
beta = M_.params( 4 );
M_.params( 5 ) = 0.35;
alpha = M_.params( 5 );
sigma = 0.007;
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 4 ) = 1;
oo_.steady_state( 2 ) = 0.5;
oo_.steady_state( 3 ) = 0.5;
oo_.steady_state( 8 ) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = sigma^2;
steady;
options_.hp_filter = 1600;
options_.order = 1;
var_list_ = char();
info = stoch_simul(var_list_);
save('Initial_Model_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('Initial_Model_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('Initial_Model_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('Initial_Model_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('Initial_Model_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('Initial_Model_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('Initial_Model_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
