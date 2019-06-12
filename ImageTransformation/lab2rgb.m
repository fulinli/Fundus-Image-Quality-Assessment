% color space conversion from l\alpha\beta to rgb
%    1.7368 <=l<= 0.0068
%    0.6428 <=a<= -0.6411
%    0.1304 <=b<= -0.1302
function [rgb] = lab2rgb(lab)
[m, n, k] = size(lab);

x = reshape(lab(:,:,1), 1, m*n);
y = reshape(lab(:,:,2), 1, m*n);
z = reshape(lab(:,:,3), 1, m*n);
xyz = [x; y; z];

M = [0.5774,0.4082,0.7071;0.5774,0.4082,-0.7071;0.5774,-0.8165,0];
lms = M*xyz;
N = [4.4679,-3.5873,0.1193;-1.2186,2.3809,-0.1624;0.0497,-0.2439,1.2045];
tmp = N*lms;

rgb = zeros(m, n, k);
rgb(:,:,1) = reshape(tmp(1,:), m, n);
rgb(:,:,2) = reshape(tmp(2,:), m, n);
rgb(:,:,3) = reshape(tmp(3,:), m, n);

s = rgb>1.0;
rgb(s) = 1.0;
s = rgb<0.0;
rgb(s) = 0.0;