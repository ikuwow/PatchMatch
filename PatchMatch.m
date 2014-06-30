% PatchMatch.m

% the Matlab code of PatchMatch algorithm

% PatchMatch returns approximate nearest neighbor field (NNF).

% targetImg: An image (usually masked by NaN. NaN is lost domain)
% sourceImg: An image from which patches are extracted, same size as targetImg.
% psz: patch size (psz x psz). Default is 9. 

% NNF: contains indices of sourceImg for each indices of targetImg
function NNF = PatchMatch(targetImg, sourceImg, psz)

% set psz to default
if (nargin<3) psz = 9; end

% grayscale images only (TODO: extend to color images)
if ndims(targetImg) > 2 targetImg = rgb2gray(targetImg); end
if ndims(sourceImg) > 2 sourceImg = rgb2gray(sourceImg); end

% targetImg = double(targetImg);
% sourceImg = double(sourceImg);

%%%%%%%%%%%%%%%%%%%%
%--  Initialize  --%
%%%%%%%%%%%%%%%%%%%%

NNF = cat(3,randi([1,size(sourceImg,1)],size(targetImg,1)),randi([1,size(sourceImg,2)],size(targetImg,2)));

%%%%%%%%%%%%%%%%%%%%%
%--  Propagation  --%
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
%--  RandomSearch  --%
%%%%%%%%%%%%%%%%%%%%%%

end
