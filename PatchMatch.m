% 
% PatchMatch.m
%
% the Matlab code of PatchMatch algorithm
% PatchMatch returns approximate nearest neighbor field (NNF).
%
% author: Ikuo Degawa <degawa@tkhm.elec.keio.ac.jp>
%
% 
% Usage: [NNF, debug] = PatchMatch(targetImg, sourceImg, psz)
% 
% Inputs: 
% - targetImg: An image (usually masked by NaN. NaN is lost domain)
% - sourceImg: An image from which patches are extracted, same size as targetImg.
% - psz: patch size (psz x psz). Default is 9. 
% 
% Outputs:
% - NNF: Nearest Neighbor Field, which contains indices of sourceImg for 
%        each corresponding indices of targetImg
% - debug: debugging information.
%

function [NNF, debug] = PatchMatch(targetImg, sourceImg, psz)

% set psz to default
if (nargin<3); psz = 9; end

% grayscale images only (TODO: extend to color images)
if ~ismatrix(targetImg); targetImg = rgb2gray(targetImg); end
if ~ismatrix(sourceImg); sourceImg = rgb2gray(sourceImg); end

targetImg = double(targetImg);
sourceImg = double(sourceImg);

%%%%%%%%%%%%%%%%%%%%
%--  Initialize  --%
%%%%%%%%%%%%%%%%%%%%

%% Parameters
max_iterations = 8;
radius = 8;
% numItr = 1;
% alpha = .5;
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
% tPatch = zeros(psz);
offsets = inf(tsz(1),tsz(2));
for ii = 1:tsz(1)
  for jj = 1:tsz(2)

    ofs = targetImg_NaN(w+ii-w:w+ii+w,w+jj-w:w+jj+w)...
          - sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);

    ofs = ofs(~isnan(ofs(:)));
    offsets(ii,jj) = sum(ofs.^2)/length(ofs);

  end
end

debug.offsets_ini = offsets;
debug.NNF_ini = NNF;


%%%%%%%%%%%%%%%%%%%%%%%%
%--  MAIN ITERATION  --%
%%%%%%%%%%%%%%%%%%%%%%%%
for iteration = 1:max_iterations

is_odd = mod(iteration,2)==1;
debug.is_odd = is_odd;

%% raster scan or reverse raster scan
if is_odd % odd
    disp([num2str(iteration),'th iteration (raster scan order) start!.']);
    ii_seq = 1:tsz(1);
    jj_seq = 1:tsz(2);
    % neighbor_dest = -1;
else % even
    disp([num2str(iteration),'th iteration (reverse raster scan order) start!.']);
    ii_seq = tsz(1):(-1):1;
    jj_seq = tsz(2):(-1):1;
    % neighbor_dest = +1;
end

fprintf('0%%----------100%%\n >'); % ten %s.

for ii = ii_seq
  for jj = jj_seq

    % TODO: if offset(ii,jj) is lower than predefined threshold, continue.

    %%%%%%%%%%%%%%%%%%%%%
    %--  Propagation  --%
    %%%%%%%%%%%%%%%%%%%%%


    %% propagate from top and left
    if is_odd %odd

            % center,      top, left
        ofs = [offsets(ii,jj), Inf, Inf ];
        if ii-1>=1; ofs(2) = offsets(ii-1,jj); end
        if jj-1>=1; ofs(3) = offsets(ii,jj-1); end
        [~,idx] = min(ofs);

        % propagate from top
        if idx==2 && NNF(ii-1,jj,1)+1+w<=ssz(1) && NNF(ii-1,jj,2)+w<=ssz(2)
            NNF(ii,jj,:) = NNF(ii-1,jj,:);
            NNF(ii,jj,1) = NNF(ii,jj,1)+1;
            tmp = targetImg_NaN(w+ii-w:w+ii+w,w+jj-w:w+jj+w)...
                  - sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
            tmp = tmp(~isnan(tmp(:)));
            offsets(ii,jj) = sum(tmp(:).^2)/length(tmp(:));

        % propagate from left
        elseif idx==3 && NNF(ii,jj-1,1)<=ssz(1) && NNF(ii,jj-1,2)+1+w<=ssz(2)
            NNF(ii,jj,:) = NNF(ii,jj-1,:);
            NNF(ii,jj,2) = NNF(ii,jj,2)+1;
            tmp = targetImg_NaN(w+ii-w:w+ii+w,w+jj-w:w+jj+w)...
                  - sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
            tmp = tmp(~isnan(tmp(:)));
            offsets(ii,jj) = sum(tmp(:).^2)/length(tmp(:));
        end

    %% propagate from bottom and right
    else %even

              % center,      bottom, right
        ofs = [offsets(ii,jj), Inf, Inf ];
        if ii+1<=tsz(1); ofs(2) = offsets(ii+1,jj); end
        if jj+1<=tsz(2); ofs(3) = offsets(ii,jj+1); end
        [~,idx] = min(ofs);

        % propagate from bottom
        if idx==2 && NNF(ii+1,jj,1)-1-w>=1 && NNF(ii+1,jj,2)-w>=1
            NNF(ii,jj,:) = NNF(ii+1,jj,:);
            NNF(ii,jj,1) = NNF(ii,jj,1)-1;
            tmp = targetImg_NaN(w+ii-w:w+ii+w,w+jj-w:w+jj+w)...
                  - sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
            tmp = tmp(~isnan(tmp(:)));
            offsets(ii,jj) = sum(tmp(:).^2)/length(tmp(:));
    
        % propagate from right
        elseif idx==3 && NNF(ii,jj+1,1)-w>=1 && NNF(ii,jj+1,2)-1-w>=1
            NNF(ii,jj,:) = NNF(ii,jj+1,:);
            NNF(ii,jj,2) = NNF(ii,jj,2)-1;
            tmp = targetImg_NaN(w+ii-w:w+ii+w,w+jj-w:w+jj+w)...
                  - sourceImg(NNF(ii,jj,1)-w:NNF(ii,jj,1)+w,NNF(ii,jj,2)-w:NNF(ii,jj,2)+w);
            tmp = tmp(~isnan(tmp(:)));
            offsets(ii,jj) = sum(tmp(:).^2)/length(tmp(:));
        end

    end    


    %%%%%%%%%%%%%%%%%%%%%%
    %--  RandomSearch  --%
    %%%%%%%%%%%%%%%%%%%%%%

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

    tmp1 = targetImg_NaN(w+ii-w:w+ii+w,w+jj-w:w+jj+w) - sourceImg(iis-w:iis+w,jjs-w:jjs+w);
    tmp2 = tmp1(~isnan(tmp1(:)));

    ofs = sum(tmp2.^2)/length(tmp2);

    if ofs < offsets(ii,jj) % if finds a more relevant patch
        NNF(ii,jj,:) = [iis;jjs];
        offsets(ii,jj) = ofs;
    end

    if mod((ii-1)*tsz(2)+jj,floor(tsz(1)*tsz(2)/10))==0
        fprintf('=')
    end

  end % jj
end % ii
fprintf('>\nDone!\n');

end % iteration

debug.offsets = offsets;

end % end of function
