function R = au_rodrigues(axis, angle, slow)

% AU_RODRIGUES  Convert axis/angle representation to rotation
%               R = AU_RODRIGUES(AXIS*ANGLE)
%               R = AU_RODRIGUES(AXIS, ANGLE)

% awf, apr07
if nargin == 0
  % unit test
  disp('Testing au_rodrigues')
  axis = [1 .2 -.3]';
  axis = axis /norm(axis);
  angle = 23 / 180 * pi;
  R = au_rodrigues(axis*angle);
  Rtrue = [0.99085454066267 0.124340615391126 0.0524088791363176;-0.0962007405070338 0.92331884094085 -0.371789907729546;-0.0946186914624554 0.363347945264319 0.926836325301361];
  au_assert_equal('R','Rtrue',1e-5,1);
  au_assert_equal('det(R)','1',1e-15,1);
  au_assert_equal('R*axis','axis',1e-15,1);
  au_assert_equal('au_rodrigues([0 0 0])','eye(3)',0,1);
  
  disp('Test timing');
  tic
  N = 10000;
  for k=1:N
    R = eye(3);
  end
  t_baseline = toc;
  fprintf('Baseline: %.1f usec\n', t_baseline/N*1e6);
  
  tic
  for k=1:N
    R = au_rodrigues(axis, angle, 1);
  end
  t1 = toc-t_baseline;
  fprintf('Slow: %.1f usec\n', t1/N*1e6);
  
  tic
  for k=1:N
    R = au_rodrigues(axis, angle);
  end
  t2 = toc-t_baseline;
  fprintf('Fast: %.1f usec\n', t2/N*1e6);
  
  if t2 > t1
    disp('FAILED: Fast is slower than slow....\n');
  end
  
  clear R
  return
end

if nargin >= 2
  w = axis*angle;
else
  w = axis;
end

if nargin == 3
  % Easy to understand, and useful for deriving fast one below as follows:
  % syms w1 w2 w3 real
  % au_ccode(au_rodrigues([w1 w2 w3]))
  theta = norm(w);
  n = w / theta;
  n_x = au_cross_matrix(n);
  R = eye(3) + n_x*sin(theta) + n_x*n_x*(1-cos(theta));
else
  w1 = w(1);
  w2 = w(2);
  w3 = w(3);
  t1 = w1*w1;
  t2 = w2*w2;
  t3 = w3*w3;
  t4 = t1+t2+t3;
  if t4 == 0
    R = eye(3);
    return
  end
  t5 = sqrt(t4);
  t6 = cos(t5);
  t7 = 1.0-t6;
  t8 = 1/t4;
  t9 = t3*t8;
  t10 = t2*t8;
  t14 = sin(t5);
  t16 = 1/t5;
  t17 = t14*w3*t16;
  t19 = t8*w1;
  t20 = t7*w2*t19;
  t23 = t14*w2*t16;
  t24 = t7*w3;
  t25 = t24*t19;
  t28 = t1*t8;
  t33 = t14*w1*t16;
  t35 = t24*t8*w2;
  R = [
    1.0+t7*(-t9-t10)        -t17+t20            t23+t25
            t17+t20  1.0+t7*(-t9-t28)          -t33+t35
           -t23+t25          t33+t35   1.0+t7*(-t10-t28)
    ];
end
