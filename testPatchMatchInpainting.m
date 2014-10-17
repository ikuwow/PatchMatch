% testPatchMatchInpainting.m
%
% the test run of PatchMatch of inpainting

clear all;
close all;

psz = 9;

SaveFolderName = datestr(now,'yymmdd-HHMMSS');
mkdir('results',SaveFolderName);

diary(fullfile('results',SaveFolderName,'log.txt'));

InputImageName = 'lena.bmp'

inImg = rgb2gray(imread(InputImageName));

mask = load('~/Documents/MATLAB/AutoShared/testimages/mask512.mat');
mask = mask.text;

[NNF, debug] = PatchMatch(inImg, [], psz, mask);
disp('PatchMatch Inpainting Done!');


diary off
