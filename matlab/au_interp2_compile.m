mex -output au_interp2_debug -g au_interp2.cxx
mex -DAU_MEX_UNCHECKED au_interp2.cxx
mex -DAU_MEX_UNCHECKED -D_OPENMP -output au_interp2_omp au_interp2.cxx
