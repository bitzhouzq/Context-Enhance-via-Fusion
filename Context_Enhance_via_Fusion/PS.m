
function F = PS(img)

R = 32;
wS = 40;
alpha = 2;

% imgEng = img.*img;
% enr = sum(imgEng(:));
% img = img/enr;

[gx, gy] = GradientSobel(img);
g = sqrt(gx.^2+gy.^2);
img = img(2:end-1, 2:end-1);

% [mm, nn] = size(img);
% m = floor(mm/wS)*wS;
% n = floor(nn/wS)*wS;
% hm = floor((mm - m)/2);
% hn = floor((nn - n)/2);
% img1 = img(hm+1:m+hm, hn+1:n+hn);
% g1 = g(hm+1:m+hm, hn+1:n+hn);
% figure, imshow(uint8(img1));

img1 = img;
g1 = g;

fun1 = @(block_struct) CSFFilterF(block_struct.data, R);
inf = blockproc(img1, [wS, wS], fun1);

fun2 = @(block_struct) sum((block_struct.data(:)).^alpha);
lambda = blockproc(g1, [wS, wS], fun2);
% fun2 = @(block_struct) CSFFilterFB(block_struct.data, R);
% lambda2 = blockproc(img1, [wS, wS], fun2);
% 
% lambda = lambda1./(lambda2+eps);
F = sum(inf(:).*lambda(:))/sum(lambda(:));
end


