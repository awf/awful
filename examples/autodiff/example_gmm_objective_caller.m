function [r, grad] = example_gmm_objective_caller(x, points)
% EXAMPLE_GMM_NLL  Return negative-log-likelihood for GMM on points

if nargin == 0
  %% Test
  d = 2;
  K = 4;
  params.log_alphas = randn(K,1);
  params.means = au_map(@(i) rand(d,1), cell(K,1));
  params.inv_cov_factors = au_map(@(i) randn(d*(d+1)/2,1), cell(K,1));
  
  points = rand(d, 10);
  
  % Flatten the parameters into a vector
  x = au_deep_vectorize(params);
  
  % Call 
  r1 = example_gmm_objective_caller(x, points)
  [r2,g] = example_gmm_objective_caller(x, points);

  au_assert_equal r1 r2 1e-15
  
  au_check_derivatives(@(x) example_gmm_objective_caller(x, points), x, g)
  
  return
end

n=size(points,2);

pts_with_conditioner = [points; ones(1,n)];

% Call example_gmm_objective_mex.cxx, 
% but see example_gmm_objective for the code that generated it
if nargout <= 1
  r = example_gmm_objective_mex(repmat(x, 1, n), pts_with_conditioner, false)';
else
  r_and_g = example_gmm_objective_mex(repmat(x, 1, n), pts_with_conditioner, true)
  r = r_and_g(1,:)'
  grad = r_and_g(2:end,:)';
end

  