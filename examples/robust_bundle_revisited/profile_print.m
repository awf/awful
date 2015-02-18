function profile_print(outprefix, mfile)


if nargin == 2
  pinfo = profile('info');
  str = profview(mfile, pinfo);
  str = htmlprune(str);
  
  
else
  fullDirname = 'c:\tmp\profout';
  pinfo = profile('info');
  issys = @(filename) strncmp(filename, matlabroot, length(matlabroot));
  for k=1:length(pinfo.FunctionTable)
    if ~issys(pinfo.FunctionTable(k).FileName)
      str = profview(k, pinfo);
      
      str=htmlprune(str);
      
      
      filename = fullfile(fullDirname,sprintf('file%d.html',k));
      fid = fopen(filename,'w');
      if fid > 0
        fprintf(fid,'%s',str);
        fclose(fid);
      else
        error(message('MATLAB:profiler:UnableToOpenFile', filename));
      end
    end
  end
end


function str = htmlprune(str)

str = regexprep(str,'<a href="matlab: profview\((\d+)\);">','<a href="file$1.html">');
% The question mark makes the .* wildcard non-greedy
str = regexprep(str,'<a href="matlab:.*?>(.*?)</a>','$1');
% Remove all the forms
str = regexprep(str,'<form.*?</form>','');

insertStr = ['<body bgcolor="#F8F8F8"><strong>' getString(message('MATLAB:profiler:StaticCopyOfReport')) '</strong><p>' ...
  '<a href="file0.html">' getString(message('MATLAB:profiler:HomeUrl')) '</a><p>'];
str = strrep(str,'<body>',insertStr);
