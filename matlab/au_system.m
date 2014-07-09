function au_system(cmd, varargin)
% AU_SYSTEM   Issue system command with matlab-separated arguments
%             The matlab system command requires you to generate a command
%             string rather than passing separate arguments
%             For example,
%             system('dir', '/w', '*.m') fails oddly
%             au_system('dir', '/w', '*.m')  -> !dir "*.m"
%             

if nargin == 0
    %% Test
    au_system dir /w *.m
    au_system('dir', '/w', 'c:\program files*')
    return
end

for k=1:length(varargin)
    s = varargin{k};
    if isempty(regexp(s, '^[a-zA-Z[]{}/:\\<>,.#''~@=-+_!£$%^*()]+$', 'once'))
        if regexp(s, '"')
            error('Cannot handle " in arg...   you should just need single quotes');
        end
        s = ['"' s '"'];
    end
    cmd = [cmd ' ' s];
end
%disp(cmd)
system(cmd);
