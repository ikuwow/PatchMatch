% testNNFDoesNotUseMissingRegion.m
% 
% Input:
% - NNF : NxMx2 matrix
% - mask : NxM matrix, 0 denotes missing region, other is 


function [bin, coord] = testNNFDoesNotUseMissingRegion(NNF,mask)

bin = true;
coord = double.empty(2,0);

for ii = 1:size(NNF,1)
    for jj = 1:size(NNF,2)
        if mask(NNF(ii,jj,1),NNF(ii,jj,2))==0
            bin = false; 
            coord(:,end+1) = [ii;jj];
        end
    end
end


end % end of function
