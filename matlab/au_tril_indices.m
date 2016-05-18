function inds = au_tril_indices(dim, diagonal)
% AU_TRIL_INDICES  Return indices of elements in lower triangle
%             INDS = AU_TRIL_INDICES(DIM [, DIAGONAL])
%             Argument DIM is the size of the matrix
%             Argument DIAGONAL is as in TRIL
%             DIAGONAL = 0 includes the main diagonal
%             DIAGONAL = -1 excludes the main diagonal


if nargin == 0
    %% Test
    au_test_begin au_tril_indices
    au_test_equal au_tril_indices(1) 1
    au_test_equal au_tril_indices(2) '[1;2;4]'
    au_test_equal au_tril_indices(3) '[1;2;3;5;6;9]'
    au_test_equal au_tril_indices(3,-1) '[2;3;6]'
    au_test_end
    
    return
end

if nargin == 1
    diagonal = 0;
end

inds = find(tril(ones(dim,dim), diagonal));
