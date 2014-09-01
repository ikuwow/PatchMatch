% testPatchMatch.m
% test run of PatchMatch.m

clear all;
close all;

inImg = rgb2gray(imread('barbara.bmp'));
srcImg = rgb2gray(imread('lena.bmp'));
% mask = ones(size(inImg,1),size(inImg,2));
% mask(100:120,200:220) = NaN; 

psz = 9;
w = (psz-1)/2;

disp('Speed Test Start');
tic
[NNF, debug] = PatchMatch(inImg, srcImg, psz);
toc

reconstImg = zeros(size(inImg));

for ii = (1+w):psz:size(inImg,1)-w
    for jj = (1+w):psz:size(inImg,2)-w
        reconstImg(ii-w:ii+w,jj-w:jj+w) = srcImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
    end
end

reconstImg = uint8(reconstImg);
PSNR(double(reconstImg),double(inImg),255)
figure(1),imshow(reconstImg);



