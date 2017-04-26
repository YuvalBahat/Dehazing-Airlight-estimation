function pairsNoiseVar = calculateInterPairsPatchesNoiseVar(DClessChosenPatches,t2DividedByt1,method)

pairsNoiseVar = reshape(nanmean((DClessChosenPatches(:,2,:,:)-bsxfun(@times,shiftdim(t2DividedByt1,-1),DClessChosenPatches(:,1,:,:))).^2),...
    partialSize(DClessChosenPatches(:,1,:,:),2:4));
switch method
    case 'divideByNormI2'
% Taking the norm of difference between the DC-less NNs (after multiplying one of them by t2/t1) divided by the number of pixels, divided by ||I_2||_2                    
        pairsNoiseVar = pairsNoiseVar./reshape(sum((DClessChosenPatches(:,2,:,:)).^2),size(pairsNoiseVar));
    case 'patchesDiffVar'%The same as above, but without dividing by the norm. Meaning it is the patches' difference (after multiplying by t2/t1) variance.
    case 'divideByAveragedNorms'
        pairsNoiseVar = pairsNoiseVar./reshape((sum((DClessChosenPatches(:,2,:,:)).^2)+...
        sum(DClessChosenPatches(:,1,:,:).^2))/2,size(pairsNoiseVar));
    case 'divideByNormOfAverage'
        pairsNoiseVar = pairsNoiseVar./reshape(nanmean(((DClessChosenPatches(:,2,:,:)+bsxfun(@times,shiftdim(t2DividedByt1,-1),...
            DClessChosenPatches(:,1,:,:)))/2).^2)*size(DClessChosenPatches,1),size(pairsNoiseVar));
    case 'divideByWeightedAveragedNorms'
% Taking the norm of difference between the DC-less NNs (after multiplying one of them by t2/t1) divided by the number of pixels,...
...divided by the average of patches norms (after multiplying one of them by t2/t1) (||I_2||_2+||I_1*t2/t1||^2)/2
        pairsNoiseVar = pairsNoiseVar./reshape((sum((DClessChosenPatches(:,2,:,:)).^2)+...
        sum(bsxfun(@times,shiftdim(t2DividedByt1,-1),DClessChosenPatches(:,1,:,:)).^2))/2,...
            size(pairsNoiseVar));
    case 'secondPatchVar'
        pairsNoiseVar = shiftdim(std(DClessChosenPatches(:,2,:,:),0,1),1);
    otherwise
        error('Not Implemented')
% % Taking the norm of difference between the DC-less NNs (after multiplying one of them by t2/t1) divided by the number of pixels, divided by our initial guess of t2^2:
%                     estimatedNeighborNoiseVar = reshape(sum((DClessChosenPatches(:,2,:,:)-...
%                         bsxfun(@times,shiftdim(t2DividedByt1,-1),DClessChosenPatches(:,1,:,:))).^2),...
%                         partialSize(DClessChosenPatches(:,1,:,:),2:4))/(parameters.patchSize^2);
%                     relevantT2Indices = convert2Dto3DIndices(patchesMapDimensions,1:3,bestMatches(2,:));
%                     estimatedNeighborNoiseVar = estimatedNeighborNoiseVar./(reshape(bsxfun(@max,initialTransmissionGuess(relevantT2Indices),sqrt(eps)),size(estimatedNeighborNoiseVar)).^2);
% %  Taking the difference between the DC-less NNs (after multiplying one of them by t2/t1) as noise, than computing its var:
%                     estimatedNeighborNoiseVar = reshape(var(DClessChosenPatches(:,2,:,:)-...
%                         bsxfun(@times,shiftdim(t2DividedByt1,-1),DClessChosenPatches(:,1,:,:))),...
%                         partialSize(DClessChosenPatches(:,1,:,:),2:4));
end