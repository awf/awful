function au_makedoc

docdir = [au_root_dir '\..\doc\'];

%%
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

fn = [docdir 'MatlabAllHelp.html'];
fprintf(2, 'Writing to [%s]\n', fn);
fd = fopen(fn, 'wt');
%fprintf(fd, '<html><title>MatlabAllHelp</title><body>\n');
for mfile = (matlab_files.m)'
  thismfile = mfile{1};
  if any(cell2mat(strfind(ok_to_ignore, thismfile)))
    continue
  end
  
  h = help(thismfile);
  if regexp(h, '^[ \t\n]*AU_')
    summary_line = regexprep(h, '\n.*', '');
    fprintf(1, '%s\n', summary_line);

    h = regexprep(h, '\n', '<br/>\n');
    fprintf(fd,'<p/><hr/>\n%s\n<pre>\n%s</pre>\n\n', thismfile, h);
  else
    fprintf(2, 'Ignoring [%s]\n', thismfile);
    fprintf('Help is:\n%s\n', h);
  end
end
fprintf(fd, '</body></html>\n');
fclose(fd);
fprintf('Created [%s]\n', fn);

%%
% Make contents.m

template = fopen([docdir 'Contents_template.txt']);
fn = [au_root_dir '\Contents.m'];
fprintf(2, 'Writing to [%s]\n', fn);
fd = fopen(fn, 'wt');
fprintf(fd, 'function AUTO_GENERATED_FROM_Contents_template_txt\n');
while ~feof(template)
  l = fgetl(template);
  if regexp(l, '^#')
    continue
  end
  
  % replace !words with help text
  [startindex, endindex] = regexp(l, '![A-Za-z0-9_]+');
  if startindex
    word = l(startindex+1:endindex);
    h = help(word);
    summary_line = regexprep(h, '\n.*', '');
    l = [l(1:startindex-1) summary_line l(endindex+1:end)];
  end
  
  fprintf(fd, '%% %s\n', l);
end  
fclose(fd);
