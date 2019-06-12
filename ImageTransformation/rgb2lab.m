% color space conversion from rgb to l\alpha\beta
%    1.7368 <=l<= 0.0068
%    0.6428 <=a<= -0.6411
%    0.1304 <=b<= -0.1302
function [lab] = rgb2lab(rgb)
% convert the RGB image into lab (l\alpha\beta) color space
[m, n, k] = size(rgb);

x = reshape(rgb(:,:,1), 1, m*n);
y = reshape(rgb(:,:,2), 1, m*n);
z = reshape(rgb(:,:,3), 1, m*n);
xyz = [x; y; z];

M = [0.3811,0.5783,0.0402;0.1967,0.7244,0.0782;0.0241,0.1288,0.8444];
lms = M*xyz;

N = [0.5774,0.5774,0.5774;0.4082,0.4082,-0.8165;0.7071,-0.7071,0];
tmp = N*lms;
lab = zeros(m, n, k);
lab(:,:,1) = reshape(tmp(1,:), m, n);
lab(:,:,2) = reshape(tmp(2,:), m, n);
lab(:,:,3) = reshape(tmp(3,:), m, n);
