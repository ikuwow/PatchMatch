% testPatchMatch.m
% test run of PatchMatch.m

clear all;
close all;

inImg = imread('barbara.bmp');
srcImg = imread('lena.bmp');
% mask = ones(size(inImg,1),size(inImg,2));
% mask(100:120,200:220) = NaN; 

psz = 9;

tic
[NNF, debug] = PatchMatch(inImg, srcImg, psz);
toc

