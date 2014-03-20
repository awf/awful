function au_makedoc

matlab_files = what(au_root_dir);

for mfile = (matlab_files.m)'
    thismfile = mfile{1};
    if ~strcmp(thismfile, 'Contents.m')
        fprintf(2,'%s\n', thismfile);
        help(thismfile)
    end
end
