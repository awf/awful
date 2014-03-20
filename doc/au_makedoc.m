function au_makedoc

matlab_files = what(au_root_dir);

for mfile = (matlab_files.m)'
    mfile = mfile{1};
    if ~strcmp(mfile, 'Contents.m')
        fprintf(2,'%s\n', mfile);
        help(mfile)
    end
end
