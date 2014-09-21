% inpaintWithPatchMatch.m
% inpainting example using patchmatch

clear all;
close all;

SaveFolderName = datestr(now,'yymmdd-HHMMSS');
mkdir('results',SaveFolderName);

diary(fullfile('results',SaveFolderName,'log.txt'));

InputImageName = 'lena.bmp'
inImg = rgb2gray(imread(InputImageName));
%{
load textmask512.mat % mask
mask(mask>0) = 1;
mask(mask==0) = NaN;
inImg = mask.*inImg;
srcImg = inImg;
%}
SourceImageName = 'barbara.bmp'
srcImg = rgb2gray(imread(SourceImageName));

% return


psz = 9;
w = (psz-1)/2;

disp('Inpainting with PatchMatch Start!');
tic
[NNF, debug] = PatchMatch(inImg, inImg, psz);
toc
disp('Done!');


fprintf('Reconstructing Output Image... ');
reconstImg = zeros(size(inImg));
for ii = (1+w):psz:size(inImg,1)-w
    for jj = (1+w):psz:size(inImg,2)-w
        reconstImg(ii-w:ii+w,jj-w:jj+w) = srcImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
    end
end
fprintf('Done!\n');

%% reconstructed image
reconstImg = uint8(reconstImg);
figure(1),imshow(reconstImg);
imwrite(reconstImg,fullfile('results',SaveFolderName,'reconstImg.bmp'),'BMP');

diary off
