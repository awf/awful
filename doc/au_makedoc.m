function au_makedoc

matlab_files = what(au_root_dir);

ok_to_ignore = {
  'au_autodiff_example_1.m'
  'au_autodiff_example_2.m'
  'au_bsx_test.m'
  'au_ccode_test.m'
  'au_levmarq_test.m'
  'au_mex.m'
  'au_mex_test.m'
  'au_mexall.m'
  'au_mfilename_test.m'
  'au_prmat_test.m'
  'au_whist_test.m'
  'au_test_test.m'
  'au_sparse_test.m'
  'Contents.m'
  };

fd = fopen('MatlabAllHelp.html', 'wt');
%fprintf(fd, '<html><title>MatlabAllHelp</title><body>\n');
for mfile = (matlab_files.m)'
  thismfile = mfile{1};
  if any(cell2mat(strfind(ok_to_ignore, thismfile)))
    continue
  end
  
  h = help(thismfile);
  if regexp(h, '^[ \t\n]*AU_')
    h = regexprep(h, '\n', '<br/>\n');
    fprintf(fd,'<p/><hr/>\n%s\n<pre>\n%s</pre>\n\n', thismfile, h);
  else
    fprintf(2, 'Ignoring [%s]\n', thismfile);
    fprintf('Help is:\n%s\n', h);
  end
end
fprintf(fd, '</body></html>\n');
fclose(fd);
disp('Created MatlabAllHelp.html');
