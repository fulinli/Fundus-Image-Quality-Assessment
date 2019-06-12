clear all; clc;

pd = 'test_eye/055/';
src_path = [pd 'src.png'];
ref_path = [pd 'ref.png'];
ini_path = [pd 'tgt.png'];
save_path = [pd 'rst.png'];

% %if needed, use state-of-the-art method (histogram matching) to produce initial image
src = imread(src_path);
ref = imread(ref_path);
tic;
ini = ColorTransferPT(src,ref);
toc;
imwrite(ini, ini_path);

disp('Starting fidelity refinement for color transfer ...');
tic;
rst = colorTransfer(pd);
toc;

disp('Starting evaluation using gamut distance and structure PSNR ...');
Is = im2double(imread(src_path));
It = im2double(imread(ref_path));
Ix = im2double(imread(ini_path));
Io = im2double(imread(save_path));

% % compute gamut distance
do = evaluate_metric(It, Io); %output1
dx = evaluate_metric(It, Ix); %output2

disp(['Gamut Distance between output image and reference image: ' num2str(do)]);
disp(['Gamut Distance between hisgram image and reference image: ' num2str(dx)]);


% % compute structure PSNR
psnr_o = compute_psnr(Io, Is); %output1
psnr_x = compute_psnr(Ix, Is); %output2

disp(['Structure PSNR between output result and source image: ' num2str(psnr_o)]);
disp(['Structure PSNR between hisgram image and source image: ' num2str(psnr_x)]);
