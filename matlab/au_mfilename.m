function mfile = au_mfilename
% AU_MFILENAME  Return filename of caller, or "[base workspace]"

if 0
    %% Test code -- use CTRL-ENTER to run
    a = au_mfilename;
    au_test_regexp(a, '[base workspace]');
    au_mfilename_test;

end

s = dbstack;
if length(s) < 2
  mfile = '[base workspace]';
else
  mfile = s(2).name;
end
