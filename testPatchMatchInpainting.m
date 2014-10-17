% testPatchMatchInpainting.m
%
% the test run of PatchMatch of inpainting

clear all;
close all;

psz = 9;
w = (psz-1)/2;

SaveFolderName = datestr(now,'yymmdd-HHMMSS');
mkdir('results',SaveFolderName);

diary(fullfile('results',SaveFolderName,'log.txt'));

InputImageName = 'lena.bmp'

inImg = rgb2gray(imread(InputImageName));

mask = load('~/Documents/MATLAB/AutoShared/testimages/mask512.mat');
mask = mask.line;
mask(mask>0) = 1;

[NNF, debug] = PatchMatch(inImg, [], psz, mask);
disp('PatchMatch Inpainting Done!');

%% NNF
fprintf('Reconstructing Output Image... ');
reconstImg = zeros(size(inImg),'uint8');
for ii = (1+w):psz:size(inImg,1)-w
    for jj = (1+w):psz:size(inImg,2)-w
        reconstImg(ii-w:ii+w,jj-w:jj+w) = inImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
    end
end
fprintf('Reconstruction Done!\n');

%% NNF_ini
fprintf('Reconstructing Output Image (NNF_ini)... ');
reconstImg = zeros(size(inImg),'uint8');
for ii = (1+w):psz:size(inImg,1)-w
    for jj = (1+w):psz:size(inImg,2)-w
        reconstImg(ii-w:ii+w,jj-w:jj+w) = inImg(debug.NNF_ini(ii,jj,1)-w:debug.NNF_ini(ii,jj,1)+w,debug.NNF_ini(ii,jj,2)-w:debug.NNF_ini(ii,jj,2)+w);
    end
end

mask_ = mask;
mask_(mask>0) = 1;
maskedImg = inImg.*mask;

figure(1),imshow(maskedImg);
figure(2),imshow(reconstImg);

diary off
