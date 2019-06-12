function [ rstRgb ] = ColorTransferPT( srcRgb, tgtRgb )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

srcRgb = double(srcRgb)/255;
tgtRgb = double(tgtRgb)/255; 

% color space transformation: rgb -> l\alpha\beta
srcLab = rgb2lab(srcRgb);
tgtLab = rgb2lab(tgtRgb);

% histogram matching
rstLab = zeros(size(srcLab));
lmts = [1.7368,  0.0068;
	    0.6428, -0.6411;
	    0.1304, -0.1302];
units = 500;
for i=1:3
    rstLab(:,:,i) = histmatching(srcLab(:,:,i), tgtLab(:,:,i), lmts(i,2), lmts(i,1), units);
end

rstRgb = lab2rgb(rstLab);

end


