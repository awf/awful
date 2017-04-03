function au_mfilename_test(caller_name)

if nargin < 1
  caller_name = 'base workspace';
end

au_test_regexp(au_mfilename, 'au_mfilename_test')
au_test_regexp(au_mfilename(-1), caller_name)
