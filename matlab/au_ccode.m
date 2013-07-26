function out = au_ccode(symobj, filename)
% AU_CCODE Generate optimized C code from symbolic expression.
%             AU_CCODE(SYMEXPR) returns a string
%             AU_CCODE(SYMEXPR, FILENAME) writes to FILENAME
%             
%          EXAMPLE:
%             syms x y real
%             au_ccode(x^3 + x^2*y)
%             au_ccode(simplify(x^3 + x^2*y))

% Author: awf@microsoft.com

if nargin == 0
    au_ccode_test
    
    return
end

if nargin < 2
    filename = [];
end

cse = feval(symengine, 'generate::optimize', symobj);
c = feval(symengine, 'generate::C', cse);
cstring = strrep(char(c), '\n', sprintf('\n'));

% Replace "t454 = " with "double t454 ="
cstring = regexprep(cstring, '\<(\w+) =', '  double $&');

% Replace "t343[r][c] =" with "t343[c * out_rows + r];
cstring = regexprep(cstring, '\<(\w+)\[(\d+)\]\[(\d+)\] =', ...
    'out_ptr[$3 * out_rows + $2] =');

if isempty(filename)
    out = cstring;
    return
end

% If there's a filename, make it a mexFunction
if ischar(filename)
    fd = fopen(filename, 'w');
else
    fd = filename;
end

[out_rows,out_cols] = size(symobj);
vars = symvar(symobj);
nvars = length(vars);
GetVars = '';
for vi = 1:nvars
    GetVars = sprintf('%s\n  double* ptr_%s = mxGetPr(prhs[%d]);', ...
        GetVars, char(vars(vi)), vi - 1);
end

body = '  /* inner loop */';
for vi = 1:nvars
    v = char(vars(vi));
    in_v = ['in_' v];
    body = sprintf('%s\n  double %s = ptr_%s[c_in*mrows + r_in];\n', body, in_v, v);
    cstring = regexprep(cstring, ['\<' v '\>'], regexptranslate('escape', in_v));
end
body = sprintf('%s\n%s', body, cstring);

% Get Template Text
tfd = fopen('mlp_ccode_template.cpp', 'r');
template = fread(tfd, inf, 'char');
fclose(tfd);

varname = inputname(1);
if isempty(varname)
    varname = '[Anonymous expression]';
end

template = char(template');
template = strrep(template, '@VarName', varname);
template = strrep(template, '@NVars', num2str(nvars));
template = strrep(template, '@Body', body);
template = strrep(template, '@OutRows', num2str(out_rows));
template = strrep(template, '@OutCols', num2str(out_cols));
template = strrep(template, '@GetVars', GetVars);

fprintf(fd, '%s', template);
if ischar(filename)
    fclose(fd);
end
