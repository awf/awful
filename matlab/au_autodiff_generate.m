function str = au_autodiff_generate(function_handle, example_arguments, example_data, filename)
% AU_AUTODIFF_GENERATE  Generate code for function and derivatives
%      To generate C code for a function and its jacobian, use
%       STR = AU_AUTODIFF_GENERATE(FUNCTION_HANDLE, EXAMPLE_ARG, EXAMPLE_DATA)
%      or 
%       AU_AUTODIFF_GENERATE(FUNCTION_HANDLE, EXAMPLE_ARG, EXAMPLE_DATA, FILE)
%      The first form returns the core code in a string, the second wraps
%      it in a MEX wrapper, and stores it in file FILE.
%      The EXAMPLE_ARG and EXAMPLE_DATA are used to determine the size, and
%      to perform finite-difference checks just in case...
%
%      FUNCTION_HANDLE must take a column vector and return a column vector.  
%      It can also take an optional argument DATA.

% awf, dec13

if nargin == 0
    %%
    f = @(x,data) [sin(x(3)/norm(x(1:2))); 14*cos((6*x(1)-6*x(3))+x(2)^2)];
    p = [1 2 3]';
    f(p)
    
    au_autodiff_generate(f, p, [], 'c:/tmp/au_autodiff_generate_test_mex.cpp');
 
    % Check it exists and is callable
    do_jacobian = true;
    out = au_autodiff_generate_test_mex([p p p], zeros(0,3), do_jacobian);
    out_val = out(1,:)
    out_jac = out(2:end,:)

    au_test_should_fail o=au_autodiff_generate_test_mex(rand(4,1),zeros(0,1),false)
    au_test_should_fail o=au_autodiff_generate_test_mex(rand(3,1),zeros(2,1),false)

    return
end

if nargin < 4
    filename = [];
end

%% Determine sizes etc
au_assert size(example_arguments,2)==1
out = function_handle(example_arguments, example_data);
au_assert size(out,2)==1

m = size(example_arguments,1);
n = size(out,1);
md = size(example_data,1);

%% Make symbolic variables and push them through.
in = sym('x', [m 1]);
sym(in, 'real');
data = sym('data', [md 1]);
sym(data, 'real');

fprintf('au_autodiff_generate: making code for f:R^%d->R^%d\n', m, n);
out_val = function_handle(in, data);
fprintf('au_autodiff_generate: computing jacobian\n');
out_jac = jacobian(out_val, in);

out = [out_val out_jac]';

% awf xxfixme: handle symbolic sparsity find(~(out_all == 0))

%% Now wrap the computation in au_autodiff_generate_template.cpp
decls = '';
for k=1:m
    decls= [decls sprintf('    double x%d = in_ptr[%d];\n', k, k-1)];
end
for k=1:md
    decls= [decls sprintf('    double data%d = data_ptr[%d];\n', k, k-1)];
end
fprintf('au_autodiff_generate: computing ccode for do_jacobian=0\n');
BodyNoJac = [decls au_ccode(out_val)];
BodyNoJac = strrep(BodyNoJac, '[0 * out_rows + ', '[');

tic
fprintf('au_autodiff_generate: computing ccode for do_jacobian=1\n');
BodyJacobian = [decls au_ccode(out)];
BodyJacobian = strrep(BodyJacobian, '[0 * out_rows + ', '[');
codegen_time = toc;

if isempty(filename)
    str = BodyJacobian;
    return
end

%% If there is a filename, make it a mexFunction

% File descriptor or name?
if ischar(filename)
    fprintf('au_autodiff_generate: writing file [%s]\n', filename);
    fd = fopen(filename, 'w');
else
    fd = filename;
    fprintf('au_autodiff_generate: writing file [fd=%d]\n', fd);
end

body0 = '  /* inner loop do_jac=0 */';
body0 = sprintf('%s\n%s', body0, BodyNoJac);
body1 = '  /* inner loop do_jac=1 */';
body1 = sprintf('%s\n%s', body1, BodyJacobian);

% Get Template Text
tfd = fopen([au_root_dir '/au_autodiff_generate_template.cpp'], 'r');
template = fread(tfd, inf, 'char');
fclose(tfd);

varname = inputname(1);
if isempty(varname)
    varname = '[Anonymous expression]';
end

template = char(template');
template = strrep(template, '@VarName', varname);
template = strrep(template, '@BodyNoJac', body0);
template = strrep(template, '@BodyJacobian', body1);
template = strrep(template, '@InRows', num2str(m));
template = strrep(template, '@OutDim', num2str(n));
template = strrep(template, '@DataRows', num2str(md));
template = strrep(template, '@OutRows', num2str(size(out,1)));

fprintf(fd, '%s', template);
if ischar(filename)
    fclose(fd);
end

%% And mex it and test it...
if ischar(filename)
    fprintf('au_autodiff_generate: mexing... [%s]\n', filename)
    mex(['-I' au_root_dir], filename);

    [~,fn,~] = fileparts(filename);
        
    % Check we're not calling an old one on the path...
    d = dir(which(fn)); 
    mex_file_age_in_seconds = (now - datenum(d.date))*24*60*60;
    au_test_assert mex_file_age_in_seconds<5
        
    fprintf('au_autodiff_generate: testing... [%s]\n', fn)
    val_f = function_handle(example_arguments, example_data);
    val_mex = feval(fn, example_arguments(:), example_data(:), true);
    
    val_f = val_f(:)';
    au_test_equal val_f val_mex(1,:) 1e-8

    J = val_mex(2:end,:)';
    f_dc = @(x) feval(fn, x, example_data(:), false)';
    timeout = max(2, codegen_time/10);
    au_check_derivatives(f_dc, example_arguments(:), J, 1e-6, 1e-6, timeout);
end

