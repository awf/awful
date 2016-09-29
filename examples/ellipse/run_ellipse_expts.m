function run_ellipse_expts(npts)
%%
% fit ellipse to n points (i.e. determine 5 parameters of ellipse) using
% various methods.

% if 0
%   global gander_results
%   %%
%   for k=1:10000
%     try
%       demo_gander
%     catch e
%       e
%     end
%     save gander_results_marg gander_results;
%   end
%   return
% end

disp('*** DEMO_GANDER ***');
% report_results clear

randn('seed', 12)

if 0
  xp2g = @(p) [1 1 1 1 pi/180] .* drwc_par2geom(p);
  xg2p = @(g) drwc_geom2par(g(:).*[1 1 1 1 180/pi]');
end


if 0
  if nargin < 1
    npts = 2.^randi([18 20]);
  end
  
  DO_LSQ_TRUST = npts < 2e5;
  DO_LSQ_LM = npts < 2e6;
  DO_HSAMPSON = npts < 2e5;
  DO_SAMPSON = npts < 2e6;
  DO_5pN_AD_HESS = npts < 1e5;
  DO_5pN_AD_NOHESS = npts < 1e4;
  DO_5pN_AD_LBFGS   = npts < 2e5;
  DO_5pN_AD_CG = false; % superslow npts < 2e5;
  DO_5pN_FDGRAD = false;
  DO_FDGRAD_CP = npts < 1e4;
  DO_FD_MARG = npts < 1e4;
  DO_ICP = 0;
  DO_SGD            = 0;
  DO_RECORD_ERRS = 0;
  DO_ITER_PLOTS = 0;
  
elseif 1
  % ICP tests
  DO_ITER_PLOTS = 1;
  
  DO_LSQ_TRUST      = 0;
  DO_LSQ_LM         = 1;
  DO_5pN_AD_HESS    = 0;
  DO_5pN_AD_NOHESS  = 0;
  DO_5pN_AD_LBFGS   = 0;
  DO_5pN_AD_CG      = 0;
  DO_5pN_FDGRAD     = 0;
  DO_HSAMPSON       = 0;~DO_ITER_PLOTS;
  DO_SAMPSON        = 0;
  DO_FDGRAD_CP      = 0;
  DO_FD_MARG        = 0;
  DO_ICP            = 1;
  DO_SGD            = 1;
  DO_RECORD_ERRS    = 1;
  if nargin < 1
    npts = 128;
  end
else
  DO_LSQ_TRUST      = 0;
  DO_LSQ_LM         = 0;
  DO_5pN_AD_HESS    = 0;
  DO_5pN_AD_NOHESS  = 0;
  DO_5pN_AD_LBFGS   = 1;
  DO_5pN_AD_CG      = 0;
  DO_5pN_FDGRAD     = 0;
  DO_HSAMPSON       = 0;
  DO_SAMPSON        = 1;
  DO_FDGRAD_CP      = 0;
  DO_FD_MARG    = 0;
  DO_ICP            = 0;
  DO_SGD            = 0;

  DO_RECORD_ERRS = 1;
  DO_ITER_PLOTS = 1;
  if nargin < 1
    npts = 128;
  end
end
fprintf('npts = %d\n', npts);


%% Make pts
noise = 0.02;

cx = 1.0;
cy = .6;
theta = 53/180*pi;
rx = .6;
ry = .8;
geom_params_gt = [cx cy rx ry theta];

[x_gt,y_gt] = ellipse_points(geom_params_gt, linspace(0,pi*2,100));

clf
if DO_ITER_PLOTS
  subplot(122)
end
hold off
plot(x_gt,y_gt,'k','linewidth',3)
hold on

t_gt = linspace(-pi/4,pi/3,npts);
[x,y] = ellipse_points(geom_params_gt, t_gt);
arclength = [0 cumsum(sqrt(diff(x).^2 + diff(y).^2))];
pts = interp1(arclength, [x' y'], linspace(0,arclength(end),npts))';
pts = pts + randn(size(pts)) * noise;
x = pts(1,:);
y = pts(2,:);
plot(x,y,'k.')

ntest = 51;
[xtest,ytest] = ellipse_points(geom_params_gt, linspace(0,pi*2,ntest));
pts_test = [xtest; ytest];

%set(draw_conic(g2p(geom_params)), 'color','w','linestyle',':');

%%
if 0
for k=1:100
  p = fitconic_gev(pts(:, randperm(npts, min(1000, npts*3/4)))');
  set(draw_conic(p), 'color', [1 1 1]/2, 'linestyle','-');
end
axis([-.5 2 -.5 1.5]*2);
end
%%
params_init = fitconic_ell(pts');
%set(draw_conic(params_init), 'color', 'r', 'linestyle',':');

record_errs_conic_handle = draw_conic(params_init);
set(record_errs_conic_handle, 'color', [1 1 1]/2, 'linewidth', 4);

%%
geom_init = p2g(params_init);
if any(geom_init([3 4]) < 0)
  error('don''t like hyperbolae');
end

% initial estimates for t are crude, but should be OK
t = linspace(0, 2*pi, 1000);
[xs,ys] = ellipse_points(geom_init, t);
tvalues = zeros(size(x));
for i=1:length(x)
  [~, tind] = min((x(i) - xs).^2 + (y(i) - ys).^2);
  tvalues(i) = t(tind);
end

if 0
  [xs,ys] = ellipse_points(geom_init, tvalues);
  plot([x; xs], [y; ys],'r-');
  axis equal
  axis([0 2 0 1.5]);
end

probspec.MaxFunEvals = 500;

tparams_init = [geom_init tvalues];
MaxFunEvals = 500;

%%
% with lsqnonlin / Levenberg Marquardt
disp('fminunc: with lsqnonlin')
options = optimset('lsqnonlin');
options.Jacobian = 'on';
options.MaxFunEvals = MaxFunEvals ;
options.MaxIter = MaxFunEvals ;
options.TolFun = 1e-8;
options.DiffMaxChange = 1e-6;
options.DiffMinChange = 1e-6;
options.DerivativeCheck = 'off';
options.Display = 'off';

lb = [];
ub = [];

f_lm = @(p) ellipse_err_lm(p, pts);
if npts < 65
  check_jacobian(f_lm, tparams_init);
end

if DO_LSQ_TRUST
  options.Algorithm = 'trust-region-reflective';
  tic;
  [tparams,~,~,~,output] = ...
    lsqnonlin(f_lm, tparams_init, lb,ub, options);
  params = g2p(tparams(1:5));
  report_results('5+n LSQ trust', params, pts, output, toc);
  
  conic_handle = draw_conic(params);
  set(conic_handle, 'color', [0 .7 0 ], 'linewidth', 4);
end

if DO_LSQ_LM
  options.Algorithm = 'levenberg-marquardt';
  record_errs_global = [];
  options.OutputFcn = @record_errs;
  
  record_errs_iter_pause = 0.1;
  
  tic;
  [tparams,~,~,~,output] = ...
    lsqnonlin(f_lm, tparams_init, lb,ub, options);
  params = g2p(tparams(1:5));
  report_results('5+n LSQ lev-marq', params, pts, output, toc);
  
  set(draw_conic(params), 'color', 'k','linestyle','--');
  drawnow
  lsq_lm_errs = record_errs_global;
  assignin('base', 'lsq_lm_errs', lsq_lm_errs);
  
end


if DO_LSQ_LM % && TEST_EARLY_STOP
  if ~evalin('base', 'exist(''test_early_allvals'', ''var'')')
    %% do this in base...
    fprintf(2, 'NOTE: Creating test_early_* in base workspace\n');
    assignin('base', 'test_early_allvals', []); 
    assignin('base', 'test_early_params_3', []); 
    assignin('base', 'test_early_params_end', []); 
  end
  vs=2:.5:8; 
  nv = length(vs);
  taus = 1+10.^-vs;
  test_early_vals = zeros(1,nv+1);
  test_early_params = zeros(nv, 6);
  tfval = lsq_lm_errs(:,2);
  testerr = lsq_lm_errs(:,3);
  for kk=1:nv
    tau = taus(kk);
    stop_iter = max(find(tfval > min(tfval)*tau));
    test_early_vals(kk) = testerr(stop_iter); 
    test_early_params(kk,:) = lsq_lm_errs(stop_iter,4:end); 
  end
  test_early_vals(end) = testerr(end);
  assignin('base', 'test_early_vals', test_early_vals);
  assignin('base', 'test_early_taus', taus);
  assignin('base', 'test_early_params', test_early_params);
  evalin('base', 'test_early_allvals =[test_early_allvals; test_early_vals ];');
  evalin('base', 'test_early_params_3 =[test_early_params_3 ; test_early_params(2,:) ];');
  evalin('base', 'test_early_params_end =[test_early_params_end ; test_early_params(end,:) ];');

   
  if 1
    subplot(121)
    hold off;
    plot([1.5 nv+.5], [1 1], 'color', [1 1 1]/2);
    hold on;
    allvals= evalin('base','test_early_allvals');
    if size(allvals,1) > 1
      boxplot(allvals(:,1:end-1)./allvals(:,end*ones(1,end-1)), 'labels', vs);
    end
    axis([1.5 nv+.5, 0 2]);
    drawnow
  end
  
end


%% with autodiff
disp('fminunc: with autodiff')

%
HessPattern = sparse(5+npts, 5+npts);
HessPattern(1:5+npts, 1:5) = 1;
HessPattern(1:5, 1:5+npts) = 1;
HessPattern = HessPattern + speye(5+npts);
au_assert_equal nnz(HessPattern) (11*npts+25)

% fminunc with Hessian
options = optimset('fminunc');
options.GradObj = 'on';
%options.HessPattern = HessPattern;
options.TolFun = 1e-8;
options.DerivativeCheck = 'off';
options.Display = 'off';

nchk = 3;
check_hessian(@(p) ellipse_err_ad(p, pts(:,1:nchk)), tparams_init(1:5+nchk))

if DO_5pN_AD_HESS
  disp('fminunc: DO_5pN_AD_HESS')
  options.Hessian = 'on';
  options.LargeScale = 'on';
  tic;
  [tparams,~,~,output] = ...
    fminunc(@(p) ellipse_err_ad(p, pts), tparams_init, options);
  params = g2p(tparams(1:5));
  report_results('5+n ad Hess', params, pts, output, toc);
  set(draw_conic(params), 'color', 'k','linestyle','-');
  drawnow
end

% Nohess
if DO_5pN_AD_NOHESS
  disp('fminunc: DO_5pN_AD_NOHESS')
  options.Hessian = 'off';
  options.LargeScale = 'off';
  tic;
  [tparams,~,~,output] = ...
    fminunc(@(p) ellipse_err_ad(p, pts), tparams_init, options);
  params = g2p(tparams(1:5));
  report_results('5+n ad noHess', params, pts, output, toc);
  set(draw_conic(params), 'color', 'y','linestyle','--');
  drawnow
end

% minimize
if DO_5pN_AD_LBFGS
  disp('fminunc: DO_5pN_LBFGS')
  tic;
  [tparams,~,output.iterations] = ...
    minimize_lbfgsb(tparams_init, @(p) ellipse_err_ad(p, pts), 2400);
  output.fval = inf;
  output.funcCount = nan;
  
  params = g2p(tparams(1:5));
  report_results('5+n ad lbfgs', params, pts, output, toc);
  drawnow
end

% minimize
if DO_5pN_AD_CG
  disp('fminunc: DO_5pN_CG')
  tic;
  [tparams,~,output.iterations] = ...
    minimize(tparams_init, @(p) ellipse_err_ad(p, pts), 2400);
  output.fval = inf;
  output.funcCount = nan;
  
  params = g2p(tparams(1:5));
  report_results('5+n ad cg', params, pts, output, toc);
  drawnow
end

if DO_5pN_FDGRAD
  % with finite-difference fminunc
  disp('fminunc: DO_5pN_FDGRAD')
  options = optimset('fminunc');
  options.GradObj = 'off';
  options.LargeScale = 'off';
  options.MaxFunEvals = MaxFunEvals;
  options.MaxIter = MaxFunEvals;
  options.TolFun = 1e-8;
  options.Display = 'off';
  
  tic;
  [tparams,fval,exitflag,output] = ...
    fminunc(@(p) err_bundle(p, pts), tparams_init, options);
  params = g2p(tparams(1:5));
  report_results('5+n fdgrad', params, pts, output, toc);
  set(draw_conic(params), 'color', 'b','linestyle','-');
  drawnow
end

%%
% with sampson
if DO_SAMPSON || DO_HSAMPSON
  disp('fminunc: with sampson')
  options = optimset('fminunc');
  options.GradObj = 'on';
  options.DerivativeCheck = 'off';
  options.MaxFunEvals = MaxFunEvals;
  options.MaxIter = MaxFunEvals;
  options.TolFun = 1e-8;
  options.TolX = 1e-8;
  options.Display = 'off';
  
  nchk = 3;
  check_hessian(@(p) ellipse_err_sampson(p, pts(:,1:nchk)), params_init)
  if DO_HSAMPSON
    disp('fminunc: DO_HSAMPSON')
    options.Hessian = 'on';
    options.LargeScale = 'on';
    record_errs_global = [];
    options.OutputFcn = @record_errs_nongeom;
    
    tic;
    [params,fval,exitflag,output] = ...
      fminunc(@(p) ellipse_err_sampson(p, pts), params_init, options);
    report_results('Hsampson', params, pts, output, toc);
    set(draw_conic(params), 'color', [1 1 0]*.7,'linestyle','-');
    drawnow
    
    hsampson_errs = record_errs_global;
    assignin('base', 'hsampson_errs', hsampson_errs);
  end
  
  if DO_SAMPSON
    disp('fminunc: DO_SAMPSON with sampson, no Hessian')
    options.Hessian = 'off';
    % options.HessPattern = sparse([ones(5,5+n); ones(n,5) eye(n)]);
    options.LargeScale = 'off';
    
    record_errs_global = [];
    options.OutputFcn = @record_errs_nongeom;
    
    tic;
    [params,fval,~,output] = ...
      fminunc(@(p) ellipse_err_sampson(p, pts), params_init, options);
    params = params_init;
    output.iterations = nan;
    output.fval = inf;
    report_results('sampson', params, pts, output, toc);
    set(draw_conic(params), 'color', [1 1 0]*1, 'linestyle','-');
    drawnow
    
    sampson_errs = record_errs_global;
    assignin('base', 'sampson_errs', sampson_errs);
    
    if 1
      record_errs_global = [];
      options.OutputFcn = @record_errs_nongeom;
      
      tic
      [params,fval_gt,~,output] = ...
        fminunc(@(p) ellipse_err_sampson(p, pts), g2p(geom_params_gt), options);
      if fval_gt < fval*0.999
        fprintf(2, 'did not converge to GT');
        report_results('sampson **** GTi', params, pts, output, toc);
      end
      
    end
  end
end

%%
if DO_FDGRAD_CP
  disp('');
  disp('fminunc: DO_FDGRAD_CP closest-point objective')
  options = optimset('fminunc');
  options.GradObj = 'off';
  options.LargeScale = 'off';
  options.MaxFunEvals = MaxFunEvals;
  options.TolFun = 1e-8;
  options.Display = 'off';
  
  tic
  [params,fval,exitflag,output] = ...
    fminunc(@(p) err_cp(p, pts), params_init, options);
  report_results('fdgrad cp', params, pts, output, toc);
  set(draw_conic(params), 'color', 'm','linestyle','-');
  drawnow
end

%%
if DO_FD_MARG
  disp('');
  disp('fminunc: DO_FD_MARG marginalized objective')
  options = optimset('fminunc');
  options.GradObj = 'off';
  options.LargeScale = 'off';
  options.MaxFunEvals = MaxFunEvals;
  options.TolFun = 1e-8;
  options.TolX = 1e-8;
  options.DiffMaxChange = 1e-6;
  options.DiffMinChange = 1e-6;
  options.Display = 'off';
  
  record_errs_global = [];
  options.OutputFcn = @record_errs;
  
  tic
  [params,fval,~,output] = ...
    fminsearch(@(p) ellipse_err_marginalized(p, noise/100, pts), geom_init, options);
  report_results('simplex marg', g2p(params), pts, output, toc);
  set(draw_conic(params), 'color', 'm','linestyle','-');
  drawnow
  
  simplex_marg_errs = record_errs_global;
  assignin('base', 'simplex_marg_errs', simplex_marg_errs);
  
  record_errs_global = [];
  options.OutputFcn = @record_errs;
  
  tic
  [params,fval,~,output] = ...
    fminunc(@(p) ellipse_err_marginalized(p, noise/100, pts), geom_init, options);
  report_results('fdgrad marg', g2p(params), pts, output, toc);
  set(draw_conic(params), 'color', 'm','linestyle','-');
  drawnow
  
  fd_marg_errs = record_errs_global;
  assignin('base', 'fd_marg_errs', fd_marg_errs);
  
  if DO_RECORD_ERRS
    clf
    tmin = 1e-2;
    q = @(x) (x/ntest);
    loglog(tmin + lsq_lm_errs(:,1), lsq_lm_errs(:,2), '.-')
    hold on
    loglog(tmin + lsq_lm_errs(:,1), q(lsq_lm_errs(:,3)), '-')
    
    loglog(tmin + sampson_errs(:,1), sampson_errs(:,2), 'r.-')
    loglog(tmin + sampson_errs(:,1), q(sampson_errs(:,3)), 'r-')
    
    lift = @(x) (x - min(x))*1e-6 + .1;
    loglog(tmin + fd_marg_errs(:,1), lift(fd_marg_errs(:,2)), 'k.-')
    loglog(tmin + fd_marg_errs(:,1), q(fd_marg_errs(:,3)), 'k-')
  end
end

%% ICP
if DO_ICP
  disp('ICP')
  
  nchk = 23;
  fchk = @(p) err_icp(p, t(1:nchk), pts(:,1:nchk));
  check_hessian(fchk, geom_init);
  
  params = geom_init;
  
  options = optimset('fminunc');
  options.GradObj = 'on';
  options.LargeScale = 'on';
  options.Hessian = 'user-supplied';
  %options.HessPattern = HessPattern;
  options.TolFun = 1e-8;
  options.DerivativeCheck = 'off';
  options.Display = 'off';
  
  output = struct('funcCount', 0, 'iterations', 0);
  f_test_err = @(p) err_cp(g2p(p), pts_test);
  errs = [0 err_cp(g2p(params), pts) f_test_err(params)];
  tic
  for iter=1:2^20
    % Find closest points to current estimate
    [~,~,tv] = au_conic_closest_point(g2p(params), pts);
    
    trange = linspace(0, 2*pi, 1000);
    [xs,ys] = ellipse_points(params, trange);
    t = zeros(size(x));
    for i=1:length(x)
      [~, tind] = min((x(i) - xs).^2 + (y(i) - ys).^2);
      t(i) = trange(tind);
    end

    
    % Minimize 5-parameter closest-point objective
    % (Should take 2 or 3 2ndorder steps)
    f_5p_cp = @(p) err_icp(p, t, pts);
    [params, fval,~,outputk] = fminunc(f_5p_cp, params, options);
    
    if 1 && (iter < 32 || 2^round(log2(iter)) == iter)
      errs(end+1,:) = [toc fval f_test_err(params)];
      tmin = .05;
      q = @(x) (x/max(x)*.2).^1;
      if DO_ITER_PLOTS,
        subplot(121)
      end
      hold off
      loglog(tmin + lsq_lm_errs(:,1), lsq_lm_errs(:,2), '.-')
      hold on
      loglog(tmin + lsq_lm_errs(:,1), q(lsq_lm_errs(:,3)), '-')
      
      if ~DO_ITER_PLOTS
        loglog(tmin + hsampson_errs(:,1), hsampson_errs(:,2), 'r.-')
        loglog(tmin + hsampson_errs(:,1), q(hsampson_errs(:,3)), 'r-')
      end
      
      loglog(tmin + errs(:,1), errs(:,2), 'k.-')
      loglog(tmin + errs(:,1), q(errs(:,3)), 'k-')
      errcp = err_cp(g2p(params), pts);
      title(sprintf('iter=%d, err = %g', iter, errcp*100));
      drawnow
    end
    output.iterations = output.iterations + 1;
    output.funcCount = output.funcCount + outputk.funcCount;
    
    if DO_ITER_PLOTS
      draw_conic(g2p(params), record_errs_conic_handle);
      drawnow
    end
  end
  assignin('base', 'icp_errs', errs);
  report_results('5 ICP Hess', g2p(params), pts, output, toc);
end

%% sgd
if DO_SGD
  %%
  disp('sgd: Stochastic gradient descent')
  
  output = struct('funcCount', 0, 'iterations', 0);
  f_test_err = @(p) err_cp(g2p(p), pts_test);
  f_geom = @(p, x) err_cp(g2p(p), x);
  
  params = p2g(params_init);
  errs = [0 f_geom(params, pts) f_test_err(params) Inf];
  tic
  batchsize = 32;
  e_best = inf;
  p_best = [];
  nu = 0;
  v = 0*params; % Momentum
  TIMEOUT = 5*60; % seconds
  for epoch=1:2^20
    % define batches
    all_batch_inds = randperm(npts);
    for kk = 1:ceil(npts/batchsize)
      i = (kk-1)*batchsize+1:min(kk*batchsize, npts);
      batch_inds = all_batch_inds(i);
      % Compute closest points and gradients
      batch_pts = pts(:,batch_inds);
      J = au_jacobian_fd(@(x) f_geom(x, batch_pts), params);
      J = full(J);
      
      % Nesterov update
      nu = .01/(1+log(epoch)/log(10));
      mu = 0.9;
      v_prev = v; 
      v = mu * v - nu * J;
      params = params - mu * v_prev + (1 + mu) * v; 
      if any(params(3:4) < 0)
        params = params + nu*J;
      end
    end
    
    e = f_geom(params, pts); % Error on full dataset
    if e < e_best
      e_best = e;
      p_best = params;
    end
    errs(end+1,:) = [toc e f_test_err(params) e_best];

    if DO_ITER_PLOTS
      subplot(122)
      draw_conic(g2p(p_best), record_errs_conic_handle);
      subplot(121)
      hold off
      loglog(errs(:,1), errs(:,2),'.')
      hold on
      loglog(errs(:,1), errs(:,4), 'r')
      axis([0 max(errs(:,1)) min(errs(:,4)) max(errs(:,4))]);
      title(sprintf('epoch %d nu %g e %g', epoch, nu, e_best));
      drawnow
    end
    
    output.iterations = output.iterations + 1;
    output.funcCount = output.funcCount + npts;
    
    if errs(end,1) > TIMEOUT
      break
    end
    
  end
  assignin('base', 'sgd_errs', errs);
  report_results('5 sgd ', params, pts, output, toc);
end


%%
try
  delete(record_errs_conic_handle)
catch e
  disp(e);
end
%report_results('dump');

%%
  function stop = record_errs(x, optimVals, state)
    stop = record_errs_nongeom(g2p(x(1:5)), optimVals, state);
  end

  function stop = record_errs_nongeom(x, optimVals, state)
    stop = false;
    switch state
      case 'iter'
        if DO_ITER_PLOTS
          draw_conic(x, record_errs_conic_handle);
          
          drawnow
          pause(record_errs_iter_pause)
        end
        
        
        if ~DO_RECORD_ERRS
          return
        end
        
        if isfield(optimVals, 'resnorm')
          fval = optimVals.resnorm;
        else
          fval = optimVals.fval;
        end
        test_err = err_cp(x, pts_test);
        record_errs_global(end+1,:) = [toc fval test_err x];
        
      otherwise
    end
  end

%%
  function report_results(alg, params, pts, output, time)
    global gander_results
    
    if nargin == 1
      switch alg
        case 'clear'
          clear gander_results
        case 'dump'
          fprintf('** results dump **\n');
          for k=1:length(gander_results)
            results = gander_results(k);
            fprintf('result: alg %17s, n=%4d, err=%6.3f, test=%6.3f, time %6.3fsec, iters %4d, funcCount %4d\n', ...
              results.method, results.n, results.err*100, results.err_test, results.time, results.iterations, results.funcCount);
          end
        otherwise
          error('unknown action');
      end
    else
      n = numel(gander_results)+1;
      gander_results(n).method = alg;
      gander_results(n).n = size(pts,2);
      gander_results(n).hash = mcen_CalcMD5(pts);
      gander_results(n).err = err_cp(params, pts);
      gander_results(n).err_test = err_cp(params, pts_test);
      gander_results(n).iterations = output.iterations;
      gander_results(n).time = time;
      gander_results(n).funcCount = output.funcCount;
      results = gander_results(n);
      fprintf('result: %-19s, n=%d, err = %.3f, time=%.3f, iters %4d, funcCount %4d\n', ...
        results.method, results.n, results.err*100, results.time, results.iterations, results.funcCount);
    end
    
  end
end

%%
function e = err_bundle(params, pts)
geom = params(1:5);
tvalues = params(6:end);
[xs,ys] = ellipse_points(geom, tvalues);
e = sum((pts(1,:) - xs).^2 + (pts(2,:) - ys).^2);
end

%%
function [e, gradient, Hessian] = err_icp(params, tvalues, pts)
% This function is defined in terms of a function ellipse_ad_mkpts,
% as follows
% Call this function f, it is called with args thus:
%  f(p, [t1 ... tn]);    % p is 5 ellipse parameters, n is number of data
% Call ellipse_ad_mkpts g, it is called for each datum as
%  g(P) = g([p ti])
% So f is defined as
%  f = sum_i g([p ti])
% Now, we have access to g's 6x1 gradient and 6x6 Hessian,
%
%   dg/dP = [dg/dp]     d2g/dP2 = [ d2g/dp2   d2g/dtdp]
%           [dg/dt]               [ d2g/dpdt  d2g/dt2 ]
%
% So f's gradient (of size 5+n) is
%     the *sum* of the p gradients,
%  df/dp = sum_i dg{i}/dp
%
% And the Hessian is sum
%  d2f/dp2 = sum_i d2g{i}/dp2

[e, gradient, Hessian] = ellipse_err_ad([params(:); tvalues(:)], pts);

gradient = gradient(1:5);
Hessian = full(Hessian(1:5, 1:5));

end

%%
function e = err_cp(params, pts)
dists2 = au_conic_closest_point(params, pts, true);
f = dists2 < 1e20;
e = sum(dists2(f)) + 100 * sum(~f);
if any(~f)
  fprintf(2, '!');
end
end
%%
function p = g2p(g)
p = drwc_geom2par(g(:).*[1 1 1 1 180/pi]');
end

function g = p2g(p)
g = [1 1 1 1 pi/180] .* drwc_par2geom(p);
end
