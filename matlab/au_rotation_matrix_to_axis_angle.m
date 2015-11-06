function w = au_rotation_matrix_to_axis_angle(R)

if nargin == 0
  au_test_begin
  
  w = rem(randn(3,1), pi);
  R = au_rodrigues(w);
  au_test_equal w au_rotation_matrix_to_axis_angle(R) 1e-10
  w = [0 0 0]';
  R = au_rodrigues(w);
  au_test_equal w au_rotation_matrix_to_axis_angle(R) 1e-10
  w = [0 0 pi]';
  R = au_rodrigues(w);
  au_test_equal w au_rotation_matrix_to_axis_angle(R) 1e-10
  au_test_end
  clear w
  return
end

M = logm(R);
w = M([6 7 2])';


% if (trace(R) > (3 - small_number))
% 
%     inverse_sinc = 1 + (1.0 / 6.0)      * R(1,2) + ...
%                        (7.0 / 360.0)    * theta_4 +
%                        (31.0 / 15120.0) * theta_6;
% 
%     rho = 0.5 * inverse_sinc * r;
% 
% end
