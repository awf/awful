function err = example_gmm_objective(params, data)

if nargin == 0
    %%
    d = 2;
    K = 3;
    params.alphas = randn(K,1);
    params.mus = au_map(@(i) randn(d,1), cell(K,1));
    params.Ls = au_map(@(i) randn(d*(d+1)/2,1), cell(K,1));

    x = au_deep_vectorize(params);
    unvec = @(x) au_deep_unvectorize(params, x);
    data = randn(d,1);
    
    f = @(x,data) example_gmm_objective(unvec(x), data);
    f(x, data)
    
    au_autodiff_generate(f, x, data, 'example_gmm_objective_mex.cxx')
    return
end

d = size(data,1);
K = size(params.alphas,1);
alphas = params.alphas;
mus = [params.mus{:}];

lower_triangle_indices = tril(ones(d,d)) ~= 0;

Ls = zeros(d, d, K, 'like', params.Ls{1});
for k=1:K
    Lparams = params.Ls{k};
    L = zeros(d,d, 'like', Ls);
    L(lower_triangle_indices) = Lparams;
    Ls(:,:,k) = L;
end

err = example_gmm(alphas, mus, Ls, data);
