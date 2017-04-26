function reshaped = reshape4outOfMemoryResistantExcecution(A,parentingDim,hardCodedSizes)
% Inputs:
% A - The array to be reshaped
% parentingDim - [1xM] vector, M=number of dimensions in reshaped. All
% dimensions that have the same parenting dim stem from this dimension in A
% hardCodedSizes - Sizes corresponding to the dimensions of reshaped, that
% are independent of the actual sizes of A. Each parenting dim can have
% maximum 1 corresponding hardCodedSizes set to nan, wich means this
% dimension size will be calculated according to the actual size of A.
if length(parentingDim)~=length(hardCodedSizes)
    error('Mismatching input sizes');
end
[uniqueDims,~,dimIndex] = unique(parentingDim);
for curDimNum = 1:length(uniqueDims)
    correspondingHardCodedSizes = hardCodedSizes(dimIndex==curDimNum);
    if sum(isnan(correspondingHardCodedSizes))>1
        error('Cannot determine deisred dimensions');
    elseif any(isnan(correspondingHardCodedSizes))
        correspondingHardCodedSizes(isnan(correspondingHardCodedSizes)) =...
            size(A,uniqueDims(curDimNum))/prod(correspondingHardCodedSizes(~isnan(correspondingHardCodedSizes)));
        hardCodedSizes(dimIndex==curDimNum) = correspondingHardCodedSizes;
    end
end
reshaped = reshape(A,hardCodedSizes);