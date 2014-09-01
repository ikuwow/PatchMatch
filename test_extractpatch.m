
clear all;
inImg = rgb2gray(imread('barbara.bmp'));
srcImg = rgb2gray(imread('lena.bmp'));
ii = 1;jj = 1;psz = 9;w=4;

NNF = cat(3,...
    randi([1+w,size(srcImg,1)-w],size(inImg)),...
    randi([1+w,size(srcImg,2)-w],size(inImg))...
);

% one by one
tic
for ii = 1:size(inImg,1)
for ii = 1:size(inImg,2)
validPixels = logical(ones(psz));
validPixels((ii-w:ii+w)<1 | (ii-w:ii+w)>size(inImg,1),:) = false;
validPixels(:,(jj-w:jj+w)<1 | (jj-w:jj+w)>size(inImg,2)) = false;
sPatch = srcImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
sPatch = sPatch(validPixels);
end
end
toc

% firstly NaN 
tic
inImg_NaN = nan(size(inImg,1)+2*w,size(inImg,2)+2*w);
inImg_NaN(1+w:size(inImg,1)+w,1+w:size(inImg,2)+w) = inImg;
for ii = 1:size(inImg,1)
for jj = 1:size(inImg,2)
sPatch = inImg_NaN(ii-w+w:ii+w+w,jj-w+w:jj+w+w);
end
end
toc
