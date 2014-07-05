% 
% PatchMatch.m
%
% the Matlab code of PatchMatch algorithm
% PatchMatch returns approximate nearest neighbor field (NNF).
% 
%
% Usage : [NNF, debug] = PatchMatch(targetImg, sourceImg, psz)
% targetImg: An image (usually masked by NaN. NaN is lost domain)
% sourceImg: An image from which patches are extracted, same size as targetImg.
% psz: patch size (psz x psz). Default is 9. 
% NNF: contains indices of sourceImg for each indices of targetImg
% debug: debugging information.

function [NNF, debug] = PatchMatch(targetImg, sourceImg, psz)

% set psz to default
if (nargin<3) psz = 9; end

% grayscale images only (TODO: extend to color images)
if ndims(targetImg) > 2 targetImg = rgb2gray(targetImg); end
if ndims(sourceImg) > 2 sourceImg = rgb2gray(sourceImg); end

targetImg = double(targetImg);
sourceImg = double(sourceImg);

%%%%%%%%%%%%%%%%%%%%
%--  Initialize  --%
%%%%%%%%%%%%%%%%%%%%

itrNum = 0;
debug = 0;

ssz = [size(sourceImg,1),size(sourceImg,2)];
tsz = [size(targetImg,1),size(targetImg,2)];
if mod(psz,2)==1
    w = (psz-1)/2;
else
    error('psz must be odd.');
end

targetImg_NaN = nan(tsz(1)+2*w,tsz(2)+2*w);
targetImg_NaN(1+w:tsz(1)+w,1+w:tsz(2)+w) = targetImg;

% NNF indices whose patches do not lap over outer range of images
NNF = cat(3,...
    randi([1+w,ssz(1)-w],tsz),...
    randi([1+w,ssz(2)-w],tsz)...
);

% initialize offsets (what a redundant code...)
% need not calcurate offset in advance? => anyway, implement!
disp('Initalizing...');
tPatch = zeros(psz);
offsets = inf(tsz(1),tsz(2));
for ii = 1:tsz(1)
    for jj = 1:tsz(2)
        % tPatch = nan(psz);

        % tPatch = targetImg( max(1,ii-w):min(ii+w,tsz(1)),max(1,jj-w):min(jj+w,tsz(2)) );
        % sPatch = sPatch(validPixels);

        tPatch = targetImg_NaN(w+ii-w:w+ii+w,w+jj-w:w+jj+w);

        sPatch = sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
            
        ofs = tPatch(:) - sPatch(:);
        ofs = ofs(~isnan(ofs(:)));
        offsets(ii,jj) = sum(ofs.^2)/length(ofs);
    end
end

debug = offsets;
% return


%%%%%%%%%%%%%%%%%%%%%
%--  Propagation  --%
%%%%%%%%%%%%%%%%%%%%%
disp('Propagating...');

% if mod(itrNum,2)==1 

%% main loop (raster scan order)
for ii = 1:tsz(1)
    for jj = 1:tsz(2)
        % disp(sprintf('ii=%d, jj=%d',ii,jj));
        % pre_ofs = offsets(ii,jj)

        ofs = [offsets(ii,jj)];
        if ii-1>=1 
            ofs = [ofs,offsets(ii-1,jj)]; 
        else
            ofs = [ofs,Inf];
        end
        if jj-1>=1
            ofs = [ofs,offsets(ii,jj-1)];
        else
            ofs = [ofs,Inf];
        end

        [~,idx] = min(ofs);
        % propagate from left
        if idx==2 && NNF(ii-1,jj,1)+1+w<=ssz(1) && NNF(ii-1,jj,2)+w<=ssz(2)
            NNF(ii,jj,:) = NNF(ii-1,jj,:);
            NNF(ii,jj,1) = NNF(ii,jj,1)+1;
            tmp = targetImg_NaN(w+ii-w:w+ii+w,w+jj-w:w+jj+w) - sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
            tmp = tmp(~isnan(tmp(:)));
            offsets(ii,jj) = sum(tmp(:).^2)/length(tmp(:));

        % propagate from above
        elseif idx==3 && NNF(ii,jj-1,1)<=ssz(1) && NNF(ii,jj-1,2)+1+w<=ssz(2)
            NNF(ii,jj,:) = NNF(ii,jj-1,:);
            NNF(ii,jj,2) = NNF(ii,jj,2)+1;
            tmp = targetImg_NaN(w+ii-w:w+ii+w,w+jj-w:w+jj+w) - sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
            tmp = tmp(~isnan(tmp(:)));
            offsets(ii,jj) = sum(tmp(:).^2)/length(tmp(:));
        end

        % it's possible to be bigger than previous offsets
        %{
        if (offsets(ii,jj) > pre_ofs)
            error('bigger offests!?');
        end
        %}

    end
end

% else 

debug.offsets = offsets;

% even iteration ( reverse raster scan order)

%  end


%%%%%%%%%%%%%%%%%%%%%%
%--  RandomSearch  --%
%%%%%%%%%%%%%%%%%%%%%%
radius = 8;
numItr = 5;
alpha = .5;
for ii = 1:tsz(1)
    imin = max(1,ii-radius);
    imax = min(tsz(1),ii+radius);
    for jj = 1:tsz(2)
        jmin = max(2,jj-radius);
        jmax = min(tsz(2),jj+radius);
        
        for itr = 1:numItr
            % rand*radius;
            ofs(itr) = 
        end

        % min
    
    end
end




disp('RandomSearch...');
radius = 8;
numItr = 1;
alpha = .5;
count = 0;
for ii = 1:tsz(1)
    for jj = 1:tsz(2)
        if jj==1
            disp(sprintf('ii=%d, jj=%d',ii,jj));
        end

        iis_min = max(1+w,NNF(ii,jj,1)-radius);
        iis_max = min(NNF(ii,jj,1)+radius,ssz(1)-w);
        jjs_min = max(1+w,NNF(ii,jj,2)-radius);
        jjs_max = min(NNF(ii,jj,2)+radius,ssz(2)-w);

        iis = floor(rand*(iis_max-iis_min+1)) + iis_min;
        jjs = floor(rand*(jjs_max-jjs_min+1)) + jjs_min;

        while iis==NNF(ii,jj,1) && jjs==NNF(ii,jj,2) % Don't allow self-matching
            iis = floor(rand*(iis_max-iis_min+1)) + iis_min;
            jjs = floor(rand*(jjs_max-jjs_min+1)) + jjs_min;
        end

        sPatch = sourceImg(iis-w:iis+w,jjs-w:jjs+w);
        tPatch = targetImg_NaN(w+ii-w:w+ii+w,w+jj-w:w+jj+w);

        tmp1 = tPatch(:) - sPatch(:);
        tmp2 = tmp1(~isnan(tmp1));
        ofs = sum(tmp2.^2)/length(tmp2);

        if ofs < offsets
            NNF(ii,jj,:) = [iis;jjs];
            offsets(ii,jj) = ofs;
            count = count + 1;
        end

    end
end

debug.offsets = offsets;

count

end % end of function
