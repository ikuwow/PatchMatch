% testPatchMatch.m
% test run of PatchMatch.m

clear all;
close all;

SaveFolderName = datestr(now,'yymmdd-HHMMSS');
mkdir('results',SaveFolderName);

diary(fullfile('results',SaveFolderName,'log.txt'));

InputImageName = 'barbara.bmp'
SourceImageName = 'lena.bmp'

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

disp('Reconstructing...');

reconstImg = zeros(size(inImg));

for ii = (1+w):psz:size(inImg,1)-w
    for jj = (1+w):psz:size(inImg,2)-w
        reconstImg(ii-w:ii+w,jj-w:jj+w) = srcImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
    end
end

reconstImg = uint8(reconstImg);
PSNRs = PSNR(double(reconstImg),double(inImg),255);
fprintf('PSNR is %.4f\n',PSNRs);
figure(1),imshow(reconstImg);

imwrite(reconstImg,fullfile('results',SaveFolderName,'reconstructed.bmp'),'BMP');

diary off

