function au_makedoc

matlab_files = what(au_root_dir);

for mfile = (matlab_files.m)'
    mfile = mfile{1};
    if ~strcmp(mfile, 'Contents.m')
        help(mfile)
    end
end
