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

ssz = [size(sourceImg,1),size(sourceImg,2)];
tsz = [size(targetImg,1),size(targetImg,2)];
if mod(psz,2)==1
    w = (psz-1)/2;
else
    error('psz must be odd.');
end


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

        tPatch = targetImg( max(1,ii-w):min(ii+w,tsz(1)),max(1,jj-w):min(jj+w,tsz(2)) );

        % [x,y] = meshgrid(ii-w:ii+w,jj-w:jj+w); % heavy...
        % existValue = (x>=1 & y>=1 & x<=tsz(1) & y<=tsz(2));
        existValue = logical(ones(psz,psz));
        existValue((ii-w:ii+w)<1,:) = false;
        existValue(:,(jj-w:jj+w)<1) = false;
        existValue((ii-w:ii+w)>tsz(1),:) = false;
        existValue(:,(jj-w:jj+w)>tsz(2)) = false;
        % debug = existValue;return
        sPatch = sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
        sPatch = sPatch(existValue);

        % debug = existValue;return
        % for itr1 = -w: w
        %     for itr2 = -w: +w
        %         if ii+itr1>=1 && jj+itr2>=1 && ii+itr1<=tsz(1) && jj+itr2<=tsz(2)
        %         % if existValue(itr1+w+1,itr2);
        %             tPatch(itr1+w+1,itr2+w+1) = targetImg(ii+itr1,jj+itr2);
        %         else 
        %             tPatch(itr1+w+1,itr2+w+1) = NaN;
        %         end
        %     end
        % end
            
        ofs = (tPatch(:) - sPatch(:)).^2;
        ofs = ofs(~isnan(ofs(:)));
        offsets(ii,jj) = sum(ofs(:))/length(ofs);
    end
end

debug = offsets;


%%%%%%%%%%%%%%%%%%%%%
%--  Propagation  --%
%%%%%%%%%%%%%%%%%%%%%
disp('Propagating...');

% if mod(itrNum,2)==1 

%% main loop (raster scan order)
for ii = 1:tsz(1)
    for jj = 1:tsz(2)
        tPatch = targetImg( max(1,ii-w):min(ii+w,tsz(1)), max(1,jj-w):min(jj+w,tsz(2)) );
        if (ii-1>=1 && jj-1 >=1)
            % compare belows
            % offsets(ii,jj);
            % offsets(ii-1,jj);
            % offsets(ii,jj-1);

            % sPatch_self = sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
            % sPatch_left = sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
            % sPatch_up = sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);

            % tPatch_self = targetImg(ii-w:ii+w,jj-w:jj+w);
            % tPatch_left = targetImg(ii-w:ii+w,jj-w-1:jj+w-1);
            % tPatch_up = targetImg(ii-w-1:ii+w-1,jj-w:jj+w);

            % selfOffset = sum((sPatch_self(:) - tPatch_self(:)).^2);
            % upOffset = sum((sPatch_up(:) - tPatch_up(:)).^2);
            % leftOffset = sum((sPatch_left(:) - tPatch_left(:)).^2);

            % switch max([upOffset,selfOffset,leftOffset]);
            % if upOffset > selfOffset


        elseif (ii -1 >= 1) 
            % compare belows
            % offsets(ii,jj);
            % offsets(ii-1,jj);
        elseif (jj -1 >= 1) 
            % compare belows
            % offsets(ii,jj);
            % offsets(ii,jj-1);
        end
        
    end
end

% else 

% even iteration ( reverse raster scan order)

%  end

%%%%%%%%%%%%%%%%%%%%%%
%--  RandomSearch  --%
%%%%%%%%%%%%%%%%%%%%%%


end % end of function
