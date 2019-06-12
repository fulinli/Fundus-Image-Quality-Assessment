function [psnr] = compute_psnr(im1, im2)

if size(im1, 3) == 3,
    im1 = rgb2ycbcr(im1);
    im1 = im1(:, :, 1);
end

if size(im2, 3) == 3,
    im2 = rgb2ycbcr(im2);
    im2 = im2(:, :, 1);
end

im1 = extr_grad_fea(im1);
im2 = extr_grad_fea(im2);

imdff = double(im1) - double(im2);
imdff = imdff(:);

rmse = sqrt(mean(imdff.^2));

psnr = 10*log10(255.0/rmse);


end

function [imFea] = extr_grad_fea( im )

[nrow, ncol] = size(im);

imFea = zeros([nrow, ncol, 2]);

% first order gradient filters
hf1 = [-1,0,1];
vf1 = [-1,0,1]';
 
imFea(:, :, 1) = conv2(im, hf1, 'same');
imFea(:, :, 2) = conv2(im, vf1, 'same');

imFea = sqrt(imFea(:, :, 1).^2 + imFea(:, :, 2).^2);
end