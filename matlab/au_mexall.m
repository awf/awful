%script au_mexall
d = cd;
newdir = mlp_dirname(mfilename('fullpath'))
cd(newdir)
mex au_sparse.cxx
mex au_ssd.cxx
mex au_whist.cxx
cd(d)
