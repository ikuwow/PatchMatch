% testPatchMatch.m
% test run of PatchMatch.m

clear all;
close all;

SaveFolderName = datestr(now,'yymmdd-HHMMSS');
mkdir('results',SaveFolderName);

diary(fullfile('results',SaveFolderName,'log.txt'));

InputImageName = 'lena.bmp'
SourceImageName = 'barbara.bmp'

inImg = rgb2gray(imread(InputImageName));
srcImg = rgb2gray(imread(SourceImageName));


% mask = ones(size(inImg,1),size(inImg,2));
% mask(100:120,200:220) = NaN; 

psz = 9;
w = (psz-1)/2;

disp('PatchMatch Start!');
tic
[NNF, debug] = PatchMatch(inImg, srcImg, psz);
toc
disp('PatchMatch Done!');

%% reconstruction %%
fprintf('Reconstructing Initial NNF Image... ');
reconstImg_ini = zeros(size(inImg));
for ii = (1+w):psz:size(inImg,1)-w
    for jj = (1+w):psz:size(inImg,2)-w
        reconstImg_ini(ii-w:ii+w,jj-w:jj+w)...
            = srcImg(debug.NNF_ini(ii,jj,1)-w:debug.NNF_ini(ii,jj,1)+w,debug.NNF_ini(ii,jj,2)-w:debug.NNF_ini(ii,jj,2)+w);
    end
end
fprintf('Done!\n');

fprintf('Reconstructing Output Image... ');
reconstImg = zeros(size(inImg));
for ii = (1+w):psz:size(inImg,1)-w
    for jj = (1+w):psz:size(inImg,2)-w
        reconstImg(ii-w:ii+w,jj-w:jj+w) = srcImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
    end
end
fprintf('Done!\n');



reconstImg = uint8(reconstImg);
reconstImg_ini = uint8(reconstImg_ini);
PSNRs = imgPSNR(reconstImg,inImg);
fprintf('PSNR is %.4f\n',PSNRs);
figure(1),imshow(reconstImg);

imwrite(inImg,fullfile('results',SaveFolderName,'inImg.bmp'),'BMP');
imwrite(srcImg,fullfile('results',SaveFolderName,'srcImg.bmp'),'BMP');
imwrite(reconstImg_ini,fullfile('results',SaveFolderName,'reconstImg_ini.bmp'),'BMP');
imwrite(reconstImg,fullfile('results',SaveFolderName,'reconstImg.bmp'),'BMP');

diary off

