function normalizedCorrelationMap = calcImageNormalizedCorrelationMap(imcols,imSize,patchSize,type,suspectedNoiseMask)
averageOverChannelsSTDs = true;

if ~exist('type','var')
    type = 'max';
end
DClessImColsNew = bsxfun(@minus,imcols,mean(imcols));

executionStruct.functionName = 'std';
executionStruct.largeArgumentIndicator = true;
executionStruct.splittingInputDimension = 2;
executionStruct.splittedOutputsIndicator = true;
executionStruct.splittingOutputDimension = 2;
[executionInfo,STDImageNew] = outOfMemoryResistantExcecution(executionStruct,imcols);
clear imcols;

shiftSize = 2;
[shiftsVectCols,shiftsVectRows] = meshgrid(-shiftSize:shiftSize:shiftSize,-shiftSize:shiftSize:shiftSize);
shiftsVect = [shiftsVectRows(:),shiftsVectCols(:)]; shiftsVect(5,:) = [];
% if ~exist('suspectedNoiseMask','var') || mean(suspectedNoiseMask(:))<0.01
%     if exist('suspectedNoiseMask','var')
%         warning('Calculating the noise level for self NC map as median over all image patches');
%     end
%     medianImVar = median(shiftdim(STDImageNew,1)).^2;
% else
%     medianImVar = median(STDImageNew(1,patchLikeCropAnImage(suspectedNoiseMask,patchSize),:).^2);
% end
    medianImVar = percentile(shiftdim(STDImageNew,1),0.3).^2;

NCmapNew = nan([size(DClessImColsNew,2),3,8]);
    executionStruct.functionName = 'calculateMeanOfProduct';
    executionStruct.largeArgumentIndicator = [true,true];
    executionStruct.splittingInputDimension = 2;
    executionStruct.splittedOutputsIndicator = true;
    executionStruct.splittingOutputDimension = 2;
for shiftNum = 1:size(NCmapNew,3)
    nonShiftedPatchesIndices = true(imSize-patchSize+1);
    shiftedPatchesIndices = nonShiftedPatchesIndices;
    if shiftsVect(shiftNum,1)>0
        nonShiftedPatchesIndices(1:shiftSize,:) = false;
        shiftedPatchesIndices(end-shiftSize+1:end,:) = false;
    elseif shiftsVect(shiftNum,1)<0
        nonShiftedPatchesIndices(end-shiftSize+1:end,:) = false;
        shiftedPatchesIndices(1:shiftSize,:) = false;
    end
    if shiftsVect(shiftNum,2)>0
        nonShiftedPatchesIndices(:,1:shiftSize) = false;
        shiftedPatchesIndices(:,end-shiftSize+1:end) = false;
    elseif shiftsVect(shiftNum,2)<0
        nonShiftedPatchesIndices(:,end-shiftSize+1:end) = false;
        shiftedPatchesIndices(:,1:shiftSize) = false;
    end
    
    [executionInfo,DCLessProduct] =...
        outOfMemoryResistantExcecution(executionStruct,DClessImColsNew(:,nonShiftedPatchesIndices(:),:),DClessImColsNew(:,shiftedPatchesIndices(:),:));
    if executionInfo.numOfRetires>0
        executionStruct.trialNum2Start = executionInfo.numOfRetires;
    end
%     DCLessProduct = mean(DClessImColsNew(:,nonShiftedPatchesIndices(:),:).*DClessImColsNew(:,shiftedPatchesIndices(:),:));
    NCmapNew(nonShiftedPatchesIndices,:,shiftNum) = shiftdim(DCLessProduct./...
        bsxfun(@max,STDImageNew(1,nonShiftedPatchesIndices,:).*STDImageNew(1,shiftedPatchesIndices,:),reshape(medianImVar,1,1,3)),1);
end
if averageOverChannelsSTDs
    normalizedCorrelationMap = abs(mean(NCmapNew,2));
else
    normalizedCorrelationMap = min(abs(NCmapNew),[],2);
end;
switch type
    case 'max'
        normalizedCorrelationMap = max(normalizedCorrelationMap,[],3);
    case 'min'
        normalizedCorrelationMap = min(normalizedCorrelationMap,[],3);
end
normalizedCorrelationMap = reshape(normalizedCorrelationMap,imSize-patchSize+1);