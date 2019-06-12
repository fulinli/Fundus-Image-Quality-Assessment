function [rst_name] = colorTransfer(pd)
%src_name = 'src.png';
%pd = 'dataset/';
t = datestr(datetime);
src_path = [pd 'src.png'];
ref_path = [pd 'ref.png'];
ini_path = [pd 'tgt.png'];
save_path = [pd 'rst.png'];
rst_name = [t 'rst.png'];

src = imread(src_path);
src = imresize(src, [240 320]);
ref = imread(ref_path);
ini = imread(ini_path);


s = double(src)/255;
f = double(ref)/255;
r = double(ini)/255;

rst = FidelityRefinement(s,f,r,[]);
%rst = imresize(rst, 2.0);

[H,w,d]=size(rst);
for i=1:H
    for j=1:w
        if (rst(i,j,1)+rst(i,j,2)+rst(i,j,3))/3 < 0.15
            rst(i,j,1)=1;
            rst(i,j,2)=1;
            rst(i,j,3)=1;
        end
    end
end

%figure;
%imshow(rst)
imwrite(rst, save_path);

end

function [ rst ] = FidelityRefinement( src, ref, ini, options)

if( (nargin < 3) || isempty(options) )
    options.k_num = 10;
    options.isMinres = true;
    options.r    = 5;
    options.eps    = 0.1^2;
end

if(~any(strcmp(fieldnames(options) , 'k_num')))
    options.k_num = 10;
end
if(~any(strcmp(fieldnames(options) , 'isMinres')))
    options.isMinres  = true;
end
if(~any(strcmp(fieldnames(options) , 'r')))
    options.r         = 5;
end
if(~any(strcmp(fieldnames(options) , 'eps')))
    options.eps      = 0.1^2;
end


am = 20;
bm = 20;
src = padarray(src, [bm,am], 'symmetric', 'both');
ini = padarray(ini, [bm,am], 'symmetric', 'both');

[h,w,c]=size(src); 

%disp('1. Clustering tone layer ...');

X = reshape(src(:),h*w,c);

kmeans_options.K                           = options.k_num;
% kmeans_options.max_ite                     = 50;
% kmeans_options.num_threads                 = 4;

%[centroids, dis, assign, nassign, quer] = yael_kmeans(X' , kmeans_options);
[idx, C] = kmeans(X, options.k_num, 'emptyaction','singleton');

%idx = assign';
%C = centroids';

%[ idx, C ] = kmeans(X,50,'emptyaction','singleton');

IDX = reshape(idx,h,w);

R = zeros(size(src));

for i = 1:h
    for j = 1:w
        for t = 1:c
            R(i,j,t) = C(IDX(i,j),t);
        end    
    end
end



% [h1,w1,c1]=size(ini); 
% 
% X1 = reshape(ini(:),h1*w1,c1);
% 
% 
% [centroids1, dis1, assign1 , nassign1 , qerr1] = yael_kmeans(X1' , kmeans_options);,
% 
% idx1 = assign1';
% C1 = centroids1';
% 
% 
% IDX1 = reshape(idx1,h1,w1);
% 
% Rr = zeros(size(ini));
% 
% for i = 1:h1
%     for j = 1:w1
%         for t = 1:c1
%             Rr(i,j,t) = C1(IDX1(i,j),t);
%         end    
%     end
% end

%disp('2. Extracting structure map ...');

S = src ./ R;

S = max(min(S,1),0);

%disp('3. Optimizing structure map ...');

S1 = zeros(size(S));

for i = 1:c
    S1(:, :, i) = guidedfilter(src(:,:,i), S(:, :, i), options.r, options.eps);
end


% L = ComputeLMatrix(src);
% 
% S1 = zeros(size(S));
% 
% lambda = 1;
% 
% for i = 1:c
%     
%     f = S(:,:,i);
%    
%     k = h*w;
%     
%     idt = spdiags(ones(k,1), 0, k, k);
%    
%     M1 = idt +  lambda*L;
%     M2 = f(:);
% 
%     %solver
%     if(options.isMinres)
%         res =minres(M1, M2);
%     else
%         res =M1\M2;   
%     end
%     
%     S1(:,:,i) =reshape(res,h,w);
%     
% end


R1 = src ./ S1;

R1 = max(min(R1,1),0);

%disp('4. Reconstructing tone layer ...');

%imwrite(R1(bm+1:end-bm, am+1:end-am, :), 'test/exam1/r0.png');

R1 = color_transfer(R1, ref);

%imwrite(R1(bm+1:end-bm, am+1:end-am, :), 'test/exam1/r1.png');

R2 = zeros(size(R));
for i=1:c
    R2(:,:,i) = gradientGuidance(R1(:,:,i), ini(:,:,i), options.isMinres);
end

%imwrite(R2(bm+1:end-bm, am+1:end-am, :), 'test/exam1/r2.png');

%disp('5. Restoring structure map ...');

O = R2 .* S1;

rst = O(bm+1:end-bm, am+1:end-am, :);

end


function [r] = gradientGuidance(s, g, isMinres)
% parameters:
%    s   - the source image (1D)
%    g   - intermediate result produced by histogram matching (1D)
%    r  - the resultant image (1D)

am = floor(size(s,1)/4);
s = padarray(s, [0,am], 'symmetric', 'both');
g = padarray(g, [0,am], 'symmetric', 'both');

[m,n] = size(s);
k = m*n;

if (~exist('isMinres','var'))
    isMinres = false;
end     
if (isempty(isMinres))
    isMinres = false;
end  
  
lambda = 1;

A = ones(k,6);
A(:,4:6) = -A(k,1);
A(:,2) = 2*A(:,2);
A(:,5) = 2*A(:,5);
d = [-m-1,-m,-m+1,m-1,m,m+1];
Dx = spdiags(A, d, k, k);

B = ones(k,6);
tmp = [1,-1,2,-2,1,-1];
for i=1:6
    B(:,i) = tmp(i)*B(:,i);
end
d = [-m-1,-m+1,-1,1,m-1,m+1];
Dy = spdiags(B, d, k, k);

tmp = Dx'*Dx+Dy'*Dy;
idt = spdiags(ones(k,1), 0, k, k);
M1 = idt + lambda*tmp;
M2 = g(:) + lambda*tmp*s(:);

if(isMinres)
    r = minres(M1, M2);%M1\M2;
else
    r = M1\M2;
end
r = reshape(r, m, n);
r = r(:, am+1:end-am, :);
end


function q = guidedfilter(I, p, r, eps)
%   GUIDEDFILTER   O(1) time implementation of guided filter.
%
%   - guidance image: I (should be a gray-scale/single channel image)
%   - filtering input image: p (should be a gray-scale/single channel image)
%   - local window radius: r
%   - regularization parameter: eps

[hei, wid] = size(I);
N = boxfilter(ones(hei, wid), r); % the size of each local patch; N=(2r+1)^2 except for boundary pixels.

mean_I = boxfilter(I, r) ./ N;
mean_p = boxfilter(p, r) ./ N;
mean_Ip = boxfilter(I.*p, r) ./ N;
cov_Ip = mean_Ip - mean_I .* mean_p; % this is the covariance of (I, p) in each local patch.

mean_II = boxfilter(I.*I, r) ./ N;
var_I = mean_II - mean_I .* mean_I;

a = cov_Ip ./ (var_I + eps); % Eqn. (5) in the paper;
b = mean_p - a .* mean_I; % Eqn. (6) in the paper;

mean_a = boxfilter(a, r) ./ N;
mean_b = boxfilter(b, r) ./ N;

q = mean_a .* I + mean_b; % Eqn. (8) in the paper;
end

function imDst = boxfilter(imSrc, r)

%   BOXFILTER   O(1) time box filtering using cumulative sum
%
%   - Definition imDst(x, y)=sum(sum(imSrc(x-r:x+r,y-r:y+r)));
%   - Running time independent of r; 
%   - Equivalent to the function: colfilt(imSrc, [2*r+1, 2*r+1], 'sliding', @sum);
%   - But much faster.

[hei, wid] = size(imSrc);
imDst = zeros(size(imSrc));

%cumulative sum over Y axis
imCum = cumsum(imSrc, 1);
%difference over Y axis
imDst(1:r+1, :) = imCum(1+r:2*r+1, :);
imDst(r+2:hei-r, :) = imCum(2*r+2:hei, :) - imCum(1:hei-2*r-1, :);
imDst(hei-r+1:hei, :) = repmat(imCum(hei, :), [r, 1]) - imCum(hei-2*r:hei-r-1, :);

%cumulative sum over X axis
imCum = cumsum(imDst, 2);
%difference over Y axis
imDst(:, 1:r+1) = imCum(:, 1+r:2*r+1);
imDst(:, r+2:wid-r) = imCum(:, 2*r+2:wid) - imCum(:, 1:wid-2*r-1);
imDst(:, wid-r+1:wid) = repmat(imCum(:, wid), [1, r]) - imCum(:, wid-2*r:wid-r-1);
end

function Io = color_transfer(Is, It)
%COLOR_TRANSFER Transfer the color style of the target image into
%   the source image. As a result, these both images share the same
%   'look and feel'.
%
%   Input
%   -------
%   Is: the source image
%   It: the target image
%
%   Output 
%   -------
%   Io: the output image
%
%   Reference 
%   -------
%   Rang M.H. Nguyen, Seon Joo Kim, Michael S. Brown
%   Illuminant Aware Gamut-Based Color Transfer
%   Computer Graphics Forum 2014
%
%   Date
%   -------
%   Nov. 24, 2014

%addpath('weighted_GE');
global pi pt Vt
[H, W, ~] = size(Is);
[Ht,Wt,~] = size(It);

%% STEP 1: White-balancing and rotating 
% Grey-Egde algorithm to estimate illuminations of the source and target
mink_norm = 5;
sigma = 2;
kappa = 10; 
[wRs, wGs, wBs, ~] = weightedGE(Is, kappa, mink_norm, sigma);
WBs = [wRs wGs wBs];
[wRt, wGt, wBt, ~] = weightedGE(It, kappa, mink_norm, sigma);
WBt = [wRt wGt wBt];
WBs = sqrt(3)*WBs/sqrt(sum(WBs.^2));
WBt = sqrt(3)*WBt/sqrt(sum(WBt.^2));
Is = reshape(Is, [], 3)';
It = reshape(It, [], 3)';
Is = diag(1./WBs) * Is;
It = diag(1./WBt) * It;

% Rotate 2 images to B axis
Is = rotate_to_Zaxis(Is, [1 1 1]);
It = rotate_to_Zaxis(It, [1 1 1]);

%% STEP 2: Matching luminance
Is = Is';
It = It';
Is(:,3) = normalizeIntensity(Is(:,3), It(:,3), H, W);

%% STEP 3: Aligning gamuts
Ms = mean(Is, 1);
Mt = mean(It, 1);
% Shift the means to the origin
Is = Is - repmat(Ms, H*W, 1);
It = It - repmat(Mt, Ht*Wt, 1);
% Compute the convex-hull
[CHi, ~] = convhull(Is, 'simplify', true);
[CHt, Vt] = convhull(It, 'simplify', true);
idi = unique(CHi(:));
idt = unique(CHt(:));
pi = Is(idi,:);
pt = It(idt,:);
% Compute the optimal matrix to align two gamuts
x0 = [0 1 1];
options = optimset('MaxIter' ,50, 'Display', 'none', 'DiffMinChange', 0.01);
x = fminunc(@myfunoptimal,x0, options);
T = [x(2)*cos(x(1)) -x(2)*sin(x(1)) 0;
     x(3)*sin(x(1))  x(3)*cos(x(1)) 0;
     0               0              1];
 
%  T = [cos(x(1))/2 -sin(x(1))/2 0;
%      sin(x(1))/2 cos(x(1))/2 0;
%      0               0              1];
% disp('The optimal matrix to align two gamuts:');
% disp(T);
% Align two gamuts
Io = T*Is';
Mt(3) = Ms(3);
Io = Io + repmat(Mt',1, H*W);

%% STEP 4: Rotate back and undo white-balancing
Io = rotate_back_from_Zaxis(Io, [1 1 1]);

Io = diag(WBt) * Io;

Io(Io < 0) = 0;
Io(Io > 1) = 1;
Io = Io';
Io = reshape(Io, H, W, 3);
end

function Io = normalizeIntensity(Is, It, ws, hs)

% create matrix
lambda = 1;
DX = createDXForward(ws,hs);
DX = DX'*DX;

DY = createDYForward(ws,hs);
DY = DY'*DY;

D = lambda*(DX + DY);
A = speye(ws*hs) + D;
clear DX DY

t = sqrt(3);
Is = Is / t;
It = It / t;
If = histeq(Is, imhist(It));


% solve the sparse matrix
b = If + D*Is;
Io = minres(A, b); %A \ b;
Io = t * Io;

end

function DX = createDXForward(N, M)
K = N*M;
DX = spdiags([-ones(K,1) -2*ones(K,1) -ones(K,1) ones(K,1) 2*ones(K,1) ones(K,1)],...
              [-N -N+1 -N+2 N N+1 N+2],K,K);
DX(1:N,:) = 0;
DX(K-N+1:end,:) = 0;

end

function DY = createDYForward(N,M)
K = N*M;
DY = spdiags([-ones(K,1) ones(K,1) -2*ones(K,1) 2*ones(K,1) -ones(K,1) ones(K,1)],...
              [-N -N+2 -1 1 N N+2],K,K);
DY(1:N,:) = 0;
DY(K-N+1:end,:) = 0;
end

function f = myfunoptimal(x)
global pi pt Vt
if (x(2) < 0.0001) && (x(2) > -0.0001)
    x(2) = 0.0001;
end
if (x(3) < 0.0001) && (x(3) > -0.0001)
    x(3) = 0.001;
end
T = [x(2)*cos(x(1)) -x(2)*sin(x(1)) 0;
     x(3)*sin(x(1))  x(3)*cos(x(1)) 0;
     0              0               1];
    
Io = pi * T';

[~, Vo] = convhull(Io);
[~, Vtotal] = convhull([Io; pt]);
f = (Vtotal - Vt) + (Vtotal - Vo);
end

function I = rotate_to_Zaxis(I, a)
b = sqrt(a(1)^2 + a(2)^2);

Txz = [a(1)/b    a(2)/b   0;
      -a(2)/b    a(1)/b   0;
       0         0        1];

c = sqrt(a(1)^2+a(2)^2+a(3)^2);
Tz =  [a(3)/c    0        -b/c;
      0          1         0;
      b/c        0        a(3)/c];

T = Tz*Txz;
I = T*I;
end

function I = rotate_back_from_Zaxis(I, a)
b = sqrt(a(1)^2 + a(2)^2);

iTxz = [a(1)/b    -a(2)/b   0;
      a(2)/b    a(1)/b   0;
       0         0        1];

c = sqrt(a(1)^2+a(2)^2+a(3)^2);
iTz =  [a(3)/c    0        b/c;
      0          1         0;
      -b/c        0        a(3)/c];

iT = iTxz*iTz;
I = iT * I;
end


function [white_R ,white_G ,white_B, output_im] = weightedGE( input_im, kappa, mink_norm, sigma, mask_cal )

if( nargin < 2 ), kappa = 1; end
if( nargin < 3 ), mink_norm = 1; end
if( nargin < 4 ), sigma = 1; end
if( nargin < 5 ), mask_cal = zeros( size( input_im, 1 ), size( input_im, 2 ) ); end 

iter = 10;     % number of iterations

tmp_ill = [1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];   % start iteration with white illuminate estimate
final_ill = tmp_ill;

input_im = double( input_im );
tmp_image = input_im;
flag = 1;
while( iter && flag )     % iteratively improve illuminant estimate
    iter = iter - 1;
    tmp_image(:, :, 1) = tmp_image(:, :, 1) ./ ( sqrt(3)*( tmp_ill(1) ) );
    tmp_image(:, :, 2) = tmp_image(:, :, 2) ./ ( sqrt(3)*( tmp_ill(2) ) );
    tmp_image(:, :, 3) = tmp_image(:, :, 3) ./ ( sqrt(3)*( tmp_ill(3) ) );
    
    [sp_var, Rw, Gw, Bw] = compute_spvar( tmp_image, sigma );
    
    mask_zeros = max( Rw, max( Gw, Bw ) ) < eps; % exclude zero gradients
    mask_pixels = ( dilation33( double( max( tmp_image, [], 3) == 255 ) ) ); % exclude saturated pixels
    mask = logical( set_border( double( ( mask_cal | mask_pixels | mask_zeros ) == 0 ), sigma+1, 0 ) );
    
    grad_im = sqrt( Rw.^2 + Gw.^2 + Bw.^2 );
    
    weight_map = ( sp_var./( grad_im ) ).^kappa;
    weight_map( weight_map > 1 ) = 1;
    
    data_Rx = power( Rw.*( weight_map ), mink_norm );
    data_Gx = power( Gw.*( weight_map ), mink_norm );
    data_Bx = power( Bw.*( weight_map ), mink_norm );
    
    tmp_ill(1) = power( sum( data_Rx( mask(:) ) ), 1/mink_norm );
    tmp_ill(2) = power( sum( data_Gx( mask(:) ) ), 1/mink_norm );
    tmp_ill(3) = power( sum( data_Bx( mask(:) ) ), 1/mink_norm );
    
    tmp_ill = tmp_ill ./ norm( tmp_ill );
    final_ill = final_ill.*tmp_ill;
    final_ill = final_ill ./ norm( final_ill );
    if ( ( acos( tmp_ill*( 1/sqrt(3)*[1 1 1]' ) )/pi*180 ) < 0.05 )  %stop iteration if chance smaller 0.05 degree (angular error)
      flag = 0;
    end
end
white_R = final_ill(1);
white_G = final_ill(2);
white_B = final_ill(3);
if ( nargout > 1 )
  output_im(:, :, 1) = input_im(:, :, 1) ./ ( sqrt(3)*( final_ill(1) ) );
  output_im(:, :, 2) = input_im(:, :, 2) ./ ( sqrt(3)*( final_ill(2) ) );
  output_im(:, :, 3) = input_im(:, :, 3) ./ ( sqrt(3)*( final_ill(3) ) );
end
end


function out = dilation33( in, it )

if ( nargin<2 )
    it = 1;
end

hh = size( in, 1 );
ll = size( in, 2 );
out = zeros( hh, ll, 3 );

while( it > 0 )
    it = it - 1;
    out(:,:,1) = [in(2:hh, :); in(hh, :)];
    out(:,:,2) = in;
    out(:,:,3) = [in(1, :); in(1:hh-1, :)];
    out2 = max( out, [], 3 );
    out(:,:,1) = [out2(:, 2:ll), out2(:, ll)];
    out(:,:,2) = out2;
    out(:,:,3) = [out2(:, 1), out2(:, 1:ll-1)];
    out = max( out, [], 3 );
    in = out;
end
end


function out = set_border( in, width, method )
% sets border to either zero method=0,or method=1 to average
if ( nargin < 3 )
  method = 1;
end

temp = ones( size( in ) );
[y x] = ndgrid( 1:size( in, 1 ), 1:size( in, 2 ) );
temp = temp.*( ( x < size( temp, 2 ) - width + 1 ) & ( x > width ) );
temp = temp.*( ( y < size( temp, 1 ) - width + 1 ) & ( y > width ) );
out  = temp.*in;
if ( method == 1 )
  out = out + ( sum( out(:) )./sum( temp(:) ) )*( ones( size( in ) ) - temp );
end
end

  function [sp_var, Rw, Gw, Bw] = compute_spvar( im, sigma )
% Compute only the specular variant. 
% The weighting scheme returned by this function is the same
% as returned by compute_qi.m but much faster. This function 
% only calculates the specular variant (giving the best results), 
% while compute_qi.m also calculates the other weighting schemes
% used in [Gijsenij et al., PAMI 2012].
%


%split color channels
R = double( im(:,:,1) );
G = double( im(:,:,2) );
B = double( im(:,:,3) );

% computation of spatial derivatives
Rx = gDer( R, sigma, 1, 0 );
Ry = gDer( R, sigma, 0, 1 );
Rw = sqrt( Rx.^2 + Ry.^2 );
 
Gx = gDer( G, sigma, 1, 0 );
Gy = gDer( G, sigma, 0, 1 );
Gw = sqrt( Gx.^2 + Gy.^2 );

Bx = gDer( B, sigma, 1, 0 );
By = gDer( B, sigma, 0, 1 );
Bw = sqrt( Bx.^2 + By.^2 );

% Opponent_der
O3_x = ( Rx + Gx + Bx ) / sqrt(3);
O3_y = ( Ry + Gy + By ) / sqrt(3);

sp_var = sqrt( O3_x.^2 + O3_y.^2 );
end

function [H] = gDer( f, sigma, iorder, jorder )

break_off_sigma = 3.;
filtersize = floor( break_off_sigma*sigma + 0.5 );

f = fill_border( f, filtersize );

x = -filtersize:1:filtersize;

Gauss = 1/( sqrt( 2*pi )*sigma )*exp( ( x.^2 )/( -2*sigma*sigma ) );

switch( iorder )
  case 0
    Gx = Gauss/sum( Gauss );
  case 1
    Gx = -( x/sigma^2 ).*Gauss;
    Gx = Gx./( sum( sum( x.*Gx ) ) );
  case 2
    Gx = ( x.^2/sigma^4 - 1/sigma^2 ).*Gauss;
    Gx = Gx - sum( Gx )/size( x, 2 );
    Gx = Gx/sum( 0.5*x.*x.*Gx );
end
H = filter2( Gx, f );

switch( jorder )
  case 0
    Gy =  Gauss/sum( Gauss );
  case 1
    Gy  =  -( x/sigma^2 ).*Gauss;
    Gy  =  Gy./( sum( sum( x.*Gy ) ) );
  case 2
    Gy = ( x.^2/sigma^4 - 1/sigma^2 ).*Gauss;
    Gy = Gy - sum( Gy )/size( x,2 );
    Gy = Gy/sum( 0.5*x.*x.*Gy );
end
H = filter2( Gy', H );

H = H( filtersize+1:size( H, 1 ) - filtersize, filtersize + 1:size( H, 2 ) - filtersize );
end

function out = fill_border( in, bw )

hh = size( in, 1 );
ww = size( in, 2 );
dd = size( in, 3 );

if( dd == 1 )
  out = zeros( hh+bw*2, ww+bw*2 );
  
  out( 1:bw, 1:bw )                       = ones( bw, bw ).*in( 1, 1 );
  out( bw+hh+1:2*bw+hh, 1:bw )            = ones( bw, bw ).*in( hh, 1 );
  out( 1:bw, bw+1+ww:2*bw+ww )            = ones( bw, bw ).*in( 1, ww );
  out( bw+hh+1:2*bw+hh, bw+1+ww:2*bw+ww ) = ones( bw, bw ).*in( hh, ww );
  out( bw+1:bw+hh, bw+1:bw+ww )           = in;
  out( 1:bw, bw+1:bw+ww )                 = ones( bw, 1 )*in( 1, : );
  out( bw+hh+1:2*bw+hh, bw+1:bw+ww )      = ones( bw, 1 )*in( hh, : );
  out( bw+1:bw+hh, 1:bw )                 = in( :, 1 )*ones( 1, bw );
  out( bw+1:bw+hh, bw+ww+1:2*bw+ww )      = in( :, ww )*ones( 1, bw );
else
  out = zeros( hh+bw*2, ww+bw*2, dd );
  for( ii = 1:dd )
    out( 1:bw, 1:bw, ii )                       = ones( bw, bw ).*in( 1, 1, ii );
    out( bw+hh+1:2*bw+hh, 1:bw, ii )            = ones( bw, bw ).*in( hh, 1, ii );
    out( 1:bw, bw+1+ww:2*bw+ww, ii )            = ones( bw, bw ).*in( 1, ww, ii );
    out( bw+hh+1:2*bw+hh, bw+1+ww:2*bw+ww, ii ) = ones( bw, bw ).*in( hh, ww, ii );
    out( bw+1:bw+hh, bw+1:bw+ww, ii )           = in( :, :, ii );
    out( 1:bw, bw+1:bw+ww, ii )                 = ones( bw, 1 )*in( 1, :, ii );
    out( bw+hh+1:2*bw+hh, bw+1:bw+ww, ii )      = ones( bw, 1 )*in( hh, :, ii );
    out( bw+1:bw+hh, 1:bw, ii )                 = in( :, 1, ii )*ones( 1, bw );
    out( bw+1:bw+hh, bw+ww+1:2*bw+ww, ii )      = in( :, ww, ii )*ones( 1, bw );
  end
end
end
