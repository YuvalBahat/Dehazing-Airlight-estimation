function [normalizedCorrelation,targetPatchesVar] = calcNormalizedCorrelation(patches,Aindices,Bindices,numOfPixelsInPatch,...
    estimatedAnoiseSTDs,estimatedBnoiseSTDs,separatelyLimitEachSTDinNC)
% The normalized correlation (Pearson's coefficient) is calculated using
% the STD of each channel and/or level (in case of a needle) separately,
% even when color channels were coupled for the purpose of finding NNs
% (meaning that all channels were divided by a sole STD value).
A = reshape(patches(:,Aindices),[size(patches,1),size(Aindices)]);
B = reshape(patches(:,Bindices),[size(patches,1),1,size(Bindices,2)]);
A = reshape(A,[numOfPixelsInPatch,size(estimatedAnoiseSTDs,3),partialSize(A,2:3)]);
B = reshape(B,numOfPixelsInPatch,size(estimatedAnoiseSTDs,3),1,size(B,3));
estimatedBnoiseSTDs = permute(estimatedBnoiseSTDs,[1,3,4,2]);
estimatedAnoiseSTDs = permute(estimatedAnoiseSTDs,[2,3,1]);
meanLessB = bsxfun(@minus,B,mean(B)); STDofB = std(meanLessB);
meanLessSmallerA = bsxfun(@minus,A,mean(A));
varOfA = var(A);
if separatelyLimitEachSTDinNC
    smallerNormalizedCorrelation = bsxfun(@rdivide,bsxfun(@times,meanLessSmallerA,meanLessB),...
        bsxfun(@times,bsxfun(@max,sqrt(varOfA),estimatedAnoiseSTDs),bsxfun(@max,STDofB,estimatedBnoiseSTDs)));
else
    smallerNormalizedCorrelation = bsxfun(@rdivide,bsxfun(@times,meanLessSmallerA,meanLessB),...
        bsxfun(@max,bsxfun(@times,sqrt(varOfA),STDofB),bsxfun(@times,estimatedAnoiseSTDs,estimatedBnoiseSTDs)));
end
smallerNormalizedCorrelation = reshape(min(mean(smallerNormalizedCorrelation),[],2),partialSize(A,3:4));
normalizedCorrelation = smallerNormalizedCorrelation;
if nargout>1
    targetPatchesVar = shiftdim(mean(varOfA,2),2);
end