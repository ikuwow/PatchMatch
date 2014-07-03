% 
% PatchMatch.m
%
% the Matlab code of PatchMatch algorithm
% PatchMatch returns approximate nearest neighbor field (NNF).

% targetImg: An image (usually masked by NaN. NaN is lost domain)
% sourceImg: An image from which patches are extracted, same size as targetImg.
% psz: patch size (psz x psz). Default is 9. 
% NNF: contains indices of sourceImg for each indices of targetImg

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

        % validPixels = logical(ones(psz,psz));
        % validPixels((ii-w:ii+w)<1 | (ii-w:ii+w)>tsz(1),:) = false;
        % validPixels(:,(jj-w:jj+w)<1 | (jj-w:jj+w)>tsz(2)) = false;
        % debug = existValue;return

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
        try

        % tPatch = targetImg( max(1,ii-w):min(ii+w,tsz(1)), max(1,jj-w):min(jj+w,tsz(2)) );
        if (ii-1>=1 && jj-1 >=1)
            % compare belows
            % offsets(ii,jj);
            % offsets(ii-1,jj);
            % offsets(ii,jj-1);

            [~, idx] = min([offsets(ii,jj) , offsets(ii-1,jj), offsets(ii,jj-1)]);
            switch idx
                case 2
                    if NNF(ii-1,jj,1)+1+w<=ssz(1) && NNF(ii-1,jj,2)+w<=ssz(2)
                    NNF(ii,jj,:) = NNF(ii-1,jj,:);
                    NNF(ii,jj,1) = NNF(ii,jj,1)+1;
                    tmp = targetImg_NaN(w+ii-w:w+ii+w,w+jj-w:w+jj+w) - sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
                    tmp = tmp(~isnan(tmp(:)));
                    offsets(ii,jj) = sum(tmp(:).^2)/length(tmp(:));

                    % offsets(ii,jj) = offsets(ii-1,jj);
                    % NNF(ii,jj,:) = NNF(ii-1,jj,:);
                    end
                case 3
                    if NNF(ii,jj-1,1)<=ssz(1) && NNF(ii,jj-1,2)+1+w<=ssz(2)
                    NNF(ii,jj,:) = NNF(ii,jj-1,:);
                    NNF(ii,jj,2) = NNF(ii,jj,2)+1;
                    tmp = targetImg_NaN(w+ii-w:w+ii+w,w+jj-w:w+jj+w) - sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
                    tmp = tmp(~isnan(tmp(:)));
                    offsets(ii,jj) = sum(tmp(:).^2)/length(tmp(:));
                    end

                    % offsets(ii,jj) = offsets(ii,jj-1);
                    % NNF(ii,jj,:) = NNF(ii,jj-1,:);
            end

        elseif (ii -1 >= 1) 
            % compare belows
            % offsets(ii,jj);
            % offsets(ii-1,jj);
            [~, idx] = min([ offsets(ii,jj) , offsets(ii-1,jj) ]);
            switch idx
                case 2 
                    if NNF(ii-1,jj,1)+1+w<=ssz(1) && NNF(ii-1,jj,2)+w<=ssz(2)
                    NNF(ii,jj,:) = NNF(ii-1,jj,:);
                    NNF(ii,jj,1) = NNF(ii,jj,1)+1;
                    tmp = targetImg_NaN(w+ii-w:w+ii+w,w+jj-w:w+jj+w) - sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
                    tmp = tmp(~isnan(tmp(:)));
                    offsets(ii,jj) = sum(tmp(:).^2)/length(tmp(:));
                    % NNF(ii,jj,:) = NNF(ii-1,jj,:)+[-1,0];
                    % tmp = (targetImg(ii-w,jj-w)-sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w)).^2
                    % offsets(ii,jj) = sum(tmp(:))/length(tmp(:));
                    % offsets(ii,jj) = offsets(ii-1,jj);
                    % NNF(ii,jj,:) = NNF(ii-1,jj,:);
                    end
            end
        elseif (jj -1 >= 1) 
            % compare belows
            % offsets(ii,jj);
            % offsets(ii,jj-1);
            [~, idx] = min([offsets(ii,jj) , offsets(ii,jj-1) ]);
            switch idx
                case 2 
                    if NNF(ii,jj-1,1)<=ssz(1) && NNF(ii,jj-1,2)+1+w<=ssz(2)
                    NNF(ii,jj,:) = NNF(ii,jj-1,:);
                    NNF(ii,jj,2) = NNF(ii,jj,2)+1;
                    tmp = targetImg_NaN(w+ii-w:w+ii+w,w+jj-w:w+jj+w) - sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
                    tmp = tmp(~isnan(tmp(:)));
                    offsets(ii,jj) = sum(tmp(:).^2)/length(tmp(:));

                    % offsets(ii,jj) = sum(tmp(:))/length(tmp(:));
                    % offsets(ii,jj) = offsets(ii,jj-1);
                    % NNF(ii,jj,:) = NNF(ii,jj-1,:);
                    end
            end
        end % endif
        % offsets(ii,jj)
        % pause(.5);

        catch err
            ii,jj
            debug = err;
            % err
            NNF(ii,jj,:)
            error('some error');
            return
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

debug = offsets;

% even iteration ( reverse raster scan order)

%  end


%%%%%%%%%%%%%%%%%%%%%%
%--  RandomSearch  --%
%%%%%%%%%%%%%%%%%%%%%%
%{
radius = 8;
numItr = 1;
alpha = .5;
for ii = 1:tsz(1)
    imin = max(1,ii-radius);
    imax = min(tsz(1),ii+radius);
    debug = [imin,imax];return
    for jj = 1:tsz(2)
        jmin = max(1,jj-radius);
        jmax = min(tsz(2),jj+radius);
        
        ofs = [];
        for itr = 1:numItr
            iii = floor(rand*(imax-imin+1)) + imin;
            jjj = floor(rand*(jmax-jmin+1)) + jmin;
            while (iii == ii && jjj == jj) % Don't allow self-matching
                iii = floor(rand*(imax-imin+1)) + imin;
                jjj = floor(rand*(jmax-jmin+1)) + jmin;
            end
            tPatch = targetImg( max(1,ii-w):min(ii+w,tsz(1)),max(1,jj-w):min(jj+w,tsz(2)) );

            existValue = logical(ones(psz,psz));
            existValue((ii-w:ii+w)<1 | (ii-w:ii+w)>tsz(1),:) = false;
            existValue(:,(jj-w:jj+w)<1 | (jj-w:jj+w)>tsz(2)) = false;

            try
                sPatch = sourceImg(iii-w:jjj+w,jjj-w:jjj+w);
            catch err
                iii
                jjj
                error('some error');
            end
            sPatch = sPatch(existValue);
            ofs = (tPatch(:) - sPatch(:)).^2;

            tmp = (tPatch(:) - sPatch(:)).^2;
            ofs(itr) = sum(tmp(:))/length(tmp(:));

            if ofs(itr) < offsets(ii,jj)
                NNF(ii,jj,:) = [iii,jjj];
                offsets(ii,jj) = ofs(itr);
            end

        end

        % min
    
    end
end
%}






end % end of function
