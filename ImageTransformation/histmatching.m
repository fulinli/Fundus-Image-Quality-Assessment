function [rst] = histmatching(src, tgt, vmin, vmax, units)
% src, tgt - 2-dimensional data (a component of an image with 3d pixels)
%
src(src<vmin) = vmin;
src(src>vmax) = vmax;
tgt(tgt<vmin) = vmin;
tgt(tgt>vmax) = vmax;

step = (vmax-vmin)/(units-1);
x = vmin:step:vmax;

hs = getHistogram(src, vmin, vmax, step);
ht = getHistogram(tgt, vmin, vmax, step);
chs = computeCumuHistogram(hs);
cht = computeCumuHistogram(ht);


rst = src;

for i=1:size(src,1)
    for j=1:size(src,2)
        tmp = src(i,j);
        posx = floor((tmp-vmin)/step) + 1;
        posy = chs(posx);
        v = 1;
        for k=2:size(cht,1)
            if posy > cht(k)
                v = k;
            else
                break;
            end
        end
        rst(i,j) = x(v);
    end
end

function [hist] = getHistogram(data2d, lo, hi, stp)
x = lo:stp:hi;
hist = zeros(size(x,2),1);
m = size(data2d,1);
n = size(data2d,2);
for i=1:m
    for j=1:n
        tmp = data2d(i,j);
        pos = floor((tmp - lo)/stp) + 1;
        hist(pos) = hist(pos) + 1;
    end
end
hist = hist / (m*n);

function [cumuhist] = computeCumuHistogram(hist)
m = size(hist,1);
cumuhist = zeros(m,1);
cumuhist(1) = hist(1);
for i=2:m
    cumuhist(i) = cumuhist(i-1) + hist(i);
end