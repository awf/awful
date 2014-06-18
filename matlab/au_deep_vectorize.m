function x = au_deep_vectorize(obj)
% AU_DEEP_VECTORIZE Flatten arbitrary structure/cell a linear vector x.
%         x = au_deep_vectorize(obj)
%         obj1 = au_deep_unvectorize(obj, x) % use obj as a template
%         au_assert_equal obj obj1

% awf, aug13

if nargin == 0
    %% Test
    a = 2;
    b.a = randn(2,3,4);
    b.b = [1 1];
    c = {randn(3,1), randn(1,1,2)};
    b.c.d = c;
    d = cell(2,2,3);
    b1 = b;
    b1.b = [2 3];
    obj = {a,{b,c;c,b1}};
    
    x1 = au_deep_vectorize(obj);
    [obj1, n] = au_deep_unvectorize(obj, x1);
    
    fprintf('x len = %d\n', length(x1));
    au_test_equal n length(x1)
    au_test_equal obj obj1
    
    
    obj = struct;
    obj(1).a = 1;
    obj(2).a = 1.2;
    obj(3).a = 2;
    x1 = au_deep_vectorize(obj);
    [obj1, n] = au_deep_unvectorize(obj, x1);
    
    fprintf('x len = %d\n', length(x1));
    au_test_equal n length(x1)
    au_test_equal obj obj1
    
    return
end

if iscell(obj)
    % Cell array
    xs = cellfun(@(x) au_deep_vectorize(x), obj, 'uniformoutput', 0);
    x = cat(1, xs{:});
elseif isstruct(obj)
    if numel(obj) > 1
        % Struct array
        x = [];
        for k=1:numel(obj)
            x = [x; au_deep_vectorize(obj(k))];
        end
    else
        % Singleton struct
        % Use fieldnames rather than structfun as the latter creates a new structure
        x = cell2mat(au_map(@(x) au_deep_vectorize(obj.(x)), fieldnames(obj)));
    end
else
    % Doesn't make sense for non-numeric types, and until I understand
    % cat(1, uint(9), 9.1), forcing it to double
    au_assert isa(obj,'double')
    x = obj(:);
end
