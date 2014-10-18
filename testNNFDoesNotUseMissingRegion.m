% testNNFDoesNotUseMissingRegion.m
% 
% Input:
% - NNF : NxMx2 matrix
% - mask : NxM matrix, 0 denotes missing region, other is 


function [bin, coord] = testNNFDoesNotUseMissingRegion(NNF,validCenters)

coord = double.empty(2,0);

for ii = 1:size(NNF,1)
    for jj = 1:size(NNF,2)
        if validCenters(NNF(ii,jj,1),NNF(ii,jj,2))==0
            coord(:,end+1) = [ii;jj];
        end
    end
end

if ~isempty(coord)
    bin = false;
else 
    bin = true;
end


end % end of function
