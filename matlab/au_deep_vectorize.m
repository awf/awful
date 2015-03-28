function [x, unvec] = au_deep_vectorize(obj, varargin)
% AU_DEEP_VECTORIZE Flatten arbitrary structure/cell a linear vector x.
%         x = au_deep_vectorize(obj)
%         obj1 = au_deep_unvectorize(obj, x) % use obj as a template
%         au_assert_equal obj obj1
%         Two different versions of the unvectorize function.
%         If you ask for a second output argument you'll get a function
%         handle which will unvectorize the argument.   But in 2015, this
%         is a perf nightmare.  If the struct is more than say 500 deep,
%         you'll just hang your machine incorrigibly.  In that case, pass a
%         second argument which is the basename of the m-file you want
%         created which will implement the correspomding unvec.  In many
%         cases this is faster and still ok.  

% awf, aug13

if nargin == 0
    %% Test
    au_test_begin
    
    a = 2;
    b.a = randn(2,3,4);
    b.b = [1 1];
    c = {randn(3,1), randn(1,1,2)};
    b.c.d = c;
    d = cell(2,2,3);
    b1 = b;
    b1.b = [2 3];
    obj1 = {a,{b,c;c,b1}};
    
    [x1, unvec1] = au_deep_vectorize(obj1);
    obj1a = unvec1(x1);
    
    % fprintf('x len = %d\n', length(x1));
    au_test_equal obj1 obj1a 0 1

    %
    [obj1b, n] = au_deep_unvectorize(obj1, x1);
    au_test_equal n length(x1)
    au_test_equal obj1 obj1b
    
    obj2 = struct;
    obj2(1).a = 1;
    obj2(2).a = 1.2;
    obj2(3).a = 2;
    x1 = au_deep_vectorize(obj2);
    [obj2a, n] = au_deep_unvectorize(obj2, x1);
    
    % fprintf('x len = %d\n', length(x1));
    au_test_equal n length(x1)
    au_test_equal obj2a obj2
    
    sx = sym('x', size(x1));
    au_test_assert isa(sx,'sym')
    obj2b = au_deep_unvectorize(obj2, sx);
    au_test_assert isa(obj2b(1).a,'sym')
    
    au_test_end
    return
end

need_unvec = nargout > 1;

sz = size(obj);

if iscell(obj)
    % Cell array
    x = [];
    unvec = '';
    for k=1:numel(obj)
      if need_unvec
        [xk, unvec_k] = au_deep_vectorize(obj{k}, varargin{:});
      else
        xk = au_deep_vectorize(obj{k}, varargin{:});
      end
      n = numel(x);
      nk = numel(xk);
      x = [x; xk];
      if need_unvec
        unvec = [unvec; 
           unvec_k(z(n+(1:nk)))];
      end
    end
    if need_unvec
      unvec = @(z) reshape(unvec(z), sz);
    end
elseif isstruct(obj)
    fields = fieldnames(obj);
    if numel(obj) > 1
        % Struct array
        if need_unvec
          [x, unvec] = au_deep_vectorize(struct2cell(obj));
          unvec = @(z) reshape(cell2struct(unvec(z), fields, 1), sz);
        else
          x = au_deep_vectorize(struct2cell(obj));
        end
    else
        % Singleton struct
        % For unvec, create a struct call of the form
        % struct('field1', unvecs{1}(x(1:4)), 'field2', unvecs{2}(x(5:8)))
        x = [];
        if need_unvec
          unvecs = cell(1,length(fields));
          unvec_str = '';
        end
        out = zeros(10000,1);
        
        n=0;
        for k = 1:length(fields)
          if need_unvec
            [xk, unvecs{k}] = au_deep_vectorize(obj.(fields{k}));
          else
            objk = obj.(fields{k});
            xk = au_deep_vectorize(objk);
          end
          nk = length(xk);
          out(n+(1:nk)) = xk;
          % x = [x; xk];
          if need_unvec
            if k > 1
              unvec_str = [unvec_str ', '];
            end
            unvec_str = [unvec_str ...
              sprintf('''%s'', unvecs{%d}(x(%d:%d))', fields{k}, k, n+1, n+nk)];
          end
          n = n + nk;
        end
        x = out(1:n);
        if need_unvec
          eval(['unvec = @(x) au_struct(' unvec_str ');']);
        end
        %x_old = cell2mat(au_map(@(x) au_deep_vectorize(obj.(x)), fieldnames(obj)));
        %au_assert_equal x x_old
    end
else
    % Doesn't make sense for non-numeric types, and until I understand
    % cat(1, uint(9), 9.1), forcing it to double.
    % BUT SLOW, so check at random.
    if rand > .99, au_assert isnumeric(obj)||islogical(obj), end
    x = obj(:);
    if need_unvec
      unvec = @(z) reshape(z, size(obj));
    end
end
