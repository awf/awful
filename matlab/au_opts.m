function opts = au_opts(varargin)
% AU_OPTS  Easy option parsing
%          OPTS = AU_OPTS('FlagA=0', 'FlagB=3', 'Algo=foo', varargin{:})
%          is all you need to remember.   The defaults are listed first,
%          the varargs override them if necessary..
%          Any value beginning with a digit is passed to str2double, 
%          any other is left as a string.

if nargin == 0
    au_opts('FlagA=0;FlagB=3', 'Algo=foo;f=1')
    au_opts('FlagA=0','FlagB=3.1', 'Algo=foo','f=1')
    return
end

for k=1:length(varargin)
    opt = varargin{k};
    while true
        [n,e] = regexp(opt, '^(?<field>\w+)=(?<val>[^;]*)', 'names', 'end');
        opt = opt(e+2:end);
        if isempty(n)
            break
        end
        val= n(1).val;
        field = n(1).field;
        if regexp(val, '^\d')
            opts.(field) = str2double(val);
        else
            opts.(field) = val;
        end
    end
    if ~isempty(opt)
        error('Bad opts fragment [%s]', opt);
    end
end
