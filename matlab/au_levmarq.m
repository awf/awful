function [x, f, log_data] = au_levmarq(x, func, opts)

% AU_LEVMARQ    Home-grown LM with line search
%               [x, J, log] = au_levmarq(x, f, opts)
%               x:     Initial estimate
%               func:  Function to be called.
%               opts:  Algorithm options (see below)
%               For some problems, matlab's lsqnonlin will be better than 
%               this, but this is useful as a "bare bones" implementation
%               that you can modify.
%         
%               OPTIONS:
%               To get a default options structure, use
%               OPTS = AU_LEVMARQ('OPTS');
%               EDIT AU_LEVMARQ  % to see options descriptions
%
%               LOG:
%               Optional third argument LOG has rows
%                 [lambda, function_value, linmin_t, funevals]

% awf, jul13

if nargin == 0
  au_levmarq_test
  return
end

if nargin == 1 && strcmp(x, 'opts')
  opts.MaxIter = 100;      % Maximum number of outer iterations
  opts.MaxFunEvals = 1000; % Maximum numbre of function calls
  opts.Display = 'iter';   % Verbosity: none, final, final+, iter
  opts.CHECK_JACOBIAN = 3; % Seconds to spend checking derivatives.
  opts.USE_LINMIN = 0;     % Use a line search?
  opts.SCHUR_SPLIT = 0;    % Use Schur decompostion. 
                           % If n = opts.SCHUR_SPLIT, and J = [A B]
                           % with cols(A) = n, then assume B'*B
                           % is fast to invert using pcg.
  opts.USE_JTJ = 1;        % Form J'*J before solving
  opts.DECOMP_LU = 0;      % 1: backslash, 0: PCG

  % Function called each time f is called.
  opts.InnerIterFcn = @(x) [];  

  % Function called before each outer iteration,
  % just before f is called again.
  % It is passed x, and may modify it before returning,
  % E.g. to re-center.
  opts.IterStartFcn = @(x) x;
  
  % Function called after each reduction in f,
  opts.PlotFcn = @(x, fval) [];

  % How to display the scalar error
  opts.DisplayErr = @(e) sum(e.^2);

  % Levenberg-Marquardt parameters.  
  opts.LAMBDA_MIN = 1e-12; 
  opts.LAMBDA_DECREASE = 2;
  opts.LAMBDA_MAX = 1e8;
  opts.LAMBDA_INCREASE_BASE = 10;
  
  x = opts;
  return
end

if nargin < 3
  opts = au_levmarq('opts');
end

switch opts.Display
  case {'none','off'}
    VERBOSE = 0;
  case 'final'
    VERBOSE = 1;
  case 'final+'
    VERBOSE = 1.5; % Single f val per iter
  case 'iter'
    VERBOSE = 2;
  otherwise
    error(['Display [' opts.Display '] not recognized']);
end

sumsq = @(x) sum(x.^2);

% Vectorize x
x = x(:);

% Call function
nparams = length(x);

%% Optionally check the Jacobian
if opts.CHECK_JACOBIAN
    timeout = opts.CHECK_JACOBIAN;
    x_fd_check = x + rand(size(x))*1e-4;
    [~,J] = func(x_fd_check ); 
    emax = au_check_derivatives(func,x_fd_check,J,1e-5, 1e-5, timeout);
    fprintf(' emax = %g, norm(x) = %g\n', full(emax), norm(x));
end

if VERBOSE >= 1
  e = func(x);
  fprintf('au_levmarq: Beginning LM params %d, residuals %d, initial err = %g\n', ...
      nparams, numel(e), opts.DisplayErr(e));

  if ~isempty(opts.PlotFcn)
      opts.PlotFcn(x, sumsq(e));
  end

  if (VERBOSE > 1) && (VERBOSE < 2)
    fprintf('au_levmarq: ');
  end
end


%% Begin the optimization
lm_lambda = 1e-6;
log_data = [];
funevals = 0;
iter = 0;
Id = speye(nparams);
while true
  % This outer loop is called for each Jacobian computation
  % We assume that computing J is an expensive operation so 
  % try to squeeze as much out of each J as possible.

  %% Recompute jacobian

  % Call user's IterStartFcn
  if ~isempty(opts.IterStartFcn)
    xnew = opts.IterStartFcn(x);
    if ~isequal(size(xnew), size(x))
        fprintf(2, 'WARNING: IterStartFcn changed size of x.\n');
    end
  end
  
  % Call func
  [e,J] = func(x);
  funevals = funevals + 1;
  iter = iter + 1;

  % Record new f
  f = sumsq(e);

  log_data = [log_data; lm_lambda f -1 funevals];

  if (VERBOSE > 1) && (VERBOSE < 2)
    fprintf(' %g', f);
  end
  
  % Estimate Jacobian norm and preconditioner for PCG
  diag_JtJ = sum(J.^2,2);
  norm_estimate = full(mean(diag_JtJ));

  %   if ~opts.DECOMP_LU
  %     Preconditioner = spdiags(diag_JtJ,0,nparams,nparams);
  %   end

  % Do a variety of line-searches using this Jacobian,
  % varying lambda each time, as well as LAMBDA_INCREASE 
  % This attempts to get to large lambda as fast as possible
  % if we are rejecting, so that the only convergence 
  % criterion is large lambda (i.e. we stop when gradient
  % descent with a tiny step cannot reduce the function)

  % Reset LAMBDA_INCREASE
  opts.LAMBDA_INCREASE = opts.LAMBDA_INCREASE_BASE;
  
  found_reduction = 0;
  while ~found_reduction
    if VERBOSE >= 2
      fprintf('au_levmarq: iter %2d, lambda=%6.2e, norm=%4.2e, ', iter, lm_lambda, norm_estimate);
      fprintf('solve ');
    end

    % Solve (J'*J + \lambda I) x = J'e
    if opts.SCHUR_SPLIT > 0
      A = J(:,1:opts.SCHUR_SPLIT);
      B = J(:,opts.SCHUR_SPLIT+1:end);
      
      AA = A'*A; AB= A'*B; BB= B'*B;
      Ia = speye(opts.SCHUR_SPLIT);
      Ib = speye(nparams - opts.SCHUR_SPLIT);
      L = BB+lm_lambda*Ib;
      %iL = inv(L);
      LAB = L\AB';
      Jte = -J'*e;
      e1 = Jte(1:opts.SCHUR_SPLIT);
      e2 = Jte(opts.SCHUR_SPLIT+1:end);
      
      iLe2 = L\e2;
      rhs = e1 - AB*iLe2;
      
      JTJMul = @(x) AA*x + lm_lambda*x - AB*(LAB*x);
      
      [dx1, ~] = pcg(JTJMul, rhs, 1e-5, 20, AA + lm_lambda*Ia);
      
      dx2 = iLe2 - LAB*dx1;
      dx = [dx1; dx2];
    elseif opts.USE_JTJ
      AugmentedJtJ = J'*J + lm_lambda*Id;
      Jte = -(J'*e);
      dx = AugmentedJtJ \ Jte;
    elseif opts.DECOMP_LU
      dx = [J; lm_lambda*Id]\[-e; zeros(nparams,1)];
    else
      % pcg
      AugmentedJtJ = J'*J + lm_lambda*Id;
      Jte = -(J'*e);
      PCG_ITERS = 20;
      [dx, ~] = pcg(AugmentedJtJ, Jte, 1e-7 * norm_estimate, PCG_ITERS); %;, Preconditioner);
      if ~all(isfinite(dx))
        error('infinite dx...');
      end
    end
    
    if opts.USE_LINMIN
        % linmin: Line search along dx
        % Define 1D function for linmin
        f1d = @(t) f1d_aux(x + t * dx, func, opts);
        % Set options for linmin
        fminbnd_options = optimset('fminbnd');
        fminbnd_options.TolX = 1e-5;  % We don't need this search to be very precise.
        fminbnd_options.MaxFunEvals = min(120, opts.MaxFunEvals - funevals);
        fminbnd_options.Display = 'off';
        
        % Execute linmin
        [t,f_test,~,fminbnd_output] = fminbnd(f1d, .2, 10, fminbnd_options);
        x_test = x + t * dx;
        e_test = func(x_test);
        funevals = funevals + fminbnd_output.funcCount;
        if VERBOSE >= 2, fprintf('linmin [t=%4.2f], f=%d, ', t, funevals); end
        % record log data for both rejections and acceptances,
        log_data = [log_data; lm_lambda, f_test, t, funevals];

    else
        x_test = x + dx;
        e_test = func(x_test);
        f_test = sumsq(e_test);
        funevals = funevals + 1;
        % record log data for both rejections and acceptances,
        log_data = [log_data; lm_lambda, f_test, 1, funevals];

    end

    f_disp = opts.DisplayErr(e_test);
    au_assert_equal numel(f_disp) 1
    if f_test < f
      if VERBOSE >= 2, drawnow; fprintf('Accept %g\n', f_disp); end
      if ~isempty(opts.PlotFcn)
          opts.PlotFcn(x_test, f_test);
      end
      
      if lm_lambda > opts.LAMBDA_MIN
        lm_lambda = lm_lambda / opts.LAMBDA_DECREASE;
      end
      f = f_test;
      break
    else
      if VERBOSE >= 2, drawnow; fprintf('**rej* %g\n', f_disp); end
      % This means linmin failed -- the whole direction is wrong, so LAMBDA_MAX should be large
      if lm_lambda > opts.LAMBDA_MAX
        endmsg = 'exceeded lambda_max';
        found_reduction = 1;
        break
      else
        lm_lambda = lm_lambda * opts.LAMBDA_INCREASE;

        % successive rejects should be exponential, so now LAMBDA_INCREASE
        % becomes 1e1, 1e2, 1e4, 1e8 etc
        opts.LAMBDA_INCREASE = opts.LAMBDA_INCREASE^2;
        % continue
      end
    end
  end

  % terminate?
  if found_reduction
    break
  end

  if funevals > opts.MaxFunEvals
    endmsg = '>MaxFunEvals';
    break
  end
  
  % use log_data to check cgce
  if size(log_data,1) > 10
    last_f = log_data(end, 2);
    mid_f = log_data(ceil(end/2), 2);
    if abs(last_f - mid_f) < 1e-8*last_f
      endmsg = 'flatlined';
      break
    end
  end
  
  % Record new x, and loop
  x = x_test;

end

if (VERBOSE >= 1) && (VERBOSE < 3)
  if VERBOSE == 1.5, fprintf('\n'); end
  if VERBOSE >= 1, fprintf('au_levmarq:'); end
  fprintf(' done [%s], rms=%g, %d iterations, %d evals\n', endmsg, rms(e), iter, funevals);
end

%%
function e = f1d_aux(x,func,opts)
opts.InnerIterFcn(x);
r = func(x);
e = sum(r(:).^2);

%%
function r = rms(e)
r = sqrt(mean(e.^2));
