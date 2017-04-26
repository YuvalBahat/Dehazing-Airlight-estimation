function indices2keep = prune2preventRepeatedPairs(queriesIndices,queriesCorrespondingLevels,pairsPriorityOrder,allowedQueryOverlap,...
    queriesMapsDimensions,patchSize,blockedShiftedPixels,numOfActuallyDesiredPairs)
if nargin<8
    numOfActuallyDesiredPairs = length(pairsPriorityOrder);
end
initialPairsPriorityVectorLength = length(pairsPriorityOrder);
pairsPriorityOrder = pairsPriorityOrder(1:min(length(pairsPriorityOrder),20*numOfActuallyDesiredPairs));
% First, rearranging the within-pair order so that patches originating from coarser levels are in the first patch. Important for the inter-pair
% shifts stage, so that I block shifts according to the bigger patch (emanating from coarser levels):
[queriesCorrespondingLevels,withinPairOrder] = sort(queriesCorrespondingLevels,'descend');
queriesIndices = queriesIndices(bsxfun(@plus,withinPairOrder,2*(0:(size(queriesIndices,2)-1))));
indices2keep = false(1,size(queriesIndices,2));
indices2keep(pairsPriorityOrder) = true; %Pairs that were a-priori eliminated were indicated by passing a nan as their patchesValues
if allowedQueryOverlap==1 || isempty(indices2keep)
    return;
end
if nargin<6
    blockedShiftedPixels = 0;
end
if size(queriesMapsDimensions,1)>1
    [~,maxDimsIndex] = max(queriesMapsDimensions(:,1));
    largestQueriesMapDimensions = queriesMapsDimensions(maxDimsIndex,:);
    downscalingFactors = (largestQueriesMapDimensions(1)+patchSize-1)./(queriesMapsDimensions(:,1)+patchSize-1);
else
    largestQueriesMapDimensions = queriesMapsDimensions;
    downscalingFactors = 1;
end
downscalingFactors(all(queriesMapsDimensions==0,2)) = 0;
maxBlockedInterPatchesPixelsShift = round2ClosestOddNum(downscalingFactors+blockedShiftedPixels,'up'); %Rounding up to the closest odd integer

maxNumOfIndices = 0;
desiredMarginsFromPatchCenter = ceil((1-allowedQueryOverlap)*patchSize*downscalingFactors-1);   desiredMarginsFromPatchCenter(all(queriesMapsDimensions==0,2)) = 0;
pairCellNumInblockedIndices = nan(1,length(indices2keep));  pairNumMapping2withinScaleValidPairNum = nan(size(pairCellNumInblockedIndices));
blockedIndices = cell(1,size(queriesMapsDimensions,1)); totalNumOfNonNanBlockedIndices = 0;
for scaleNum = 1:size(queriesMapsDimensions,1)
    bestMatchesFinerScaleIndicator = min(queriesCorrespondingLevels)==scaleNum;
    pairCellNumInblockedIndices(bestMatchesFinerScaleIndicator & indices2keep) = scaleNum;
    pairNumMapping2withinScaleValidPairNum(bestMatchesFinerScaleIndicator & indices2keep) = 1:sum(indices2keep(bestMatchesFinerScaleIndicator));
%     queriesMapDimensions = queriesMapsDimensions(scaleNum,:);
    highestScalePatchesIndecis = rescalePatchMapLocations2HighestLevel(queriesIndices(:,bestMatchesFinerScaleIndicator),...
        queriesCorrespondingLevels(:,bestMatchesFinerScaleIndicator),queriesMapsDimensions,patchSize);
    [firstPatchRow,firstPatchCol] = ind2sub(queriesMapsDimensions(1,:),highestScalePatchesIndecis(1,:));
    [secondPatchRow,secondPatchCol] = ind2sub(queriesMapsDimensions(1,:),highestScalePatchesIndecis(2,:));
%     [firstPatchRow,firstPatchCol] = ind2sub(queriesMapDimensions,queriesIndices(1,bestMatchesFinerScaleIndicator));
%     [secondPatchRow,secondPatchCol] = ind2sub(queriesMapDimensions,queriesIndices(2,bestMatchesFinerScaleIndicator));
%     if size(queriesMapsDimensions,1)>1 && downscalingFactors(scaleNum)~=1 % Transforming row and column's patch indices to highest level size image row & col patch indices:
%         firstPatchRow = round((firstPatchRow+floor(patchSize/2))*downscalingFactors(scaleNum)-floor(patchSize/2));
%         firstPatchCol = round((firstPatchCol+floor(patchSize/2))*downscalingFactors(scaleNum)-floor(patchSize/2));
%         secondPatchRow = round((secondPatchRow+floor(patchSize/2))*downscalingFactors(scaleNum)-floor(patchSize/2));
%         secondPatchCol = round((secondPatchCol+floor(patchSize/2))*downscalingFactors(scaleNum)-floor(patchSize/2));
%     end
    % Allowing shifts in one patch of the pair
%     marginsVector = (-desiredMarginsFromPatchCenter(scaleNum):desiredMarginsFromPatchCenter(scaleNum));
%     marginsTwoDmask = bsxfun(@plus,abs(marginsVector.'),abs(marginsVector));
%     marginsTwoDmask(marginsTwoDmask>desiredMarginsFromPatchCenter(scaleNum)) = nan;
%     marginsTwoDmask(~isnan(marginsTwoDmask)) = 1;
    curScaleIndices2keep = indices2keep(bestMatchesFinerScaleIndicator);
    curScaleCorrespondingScales = queriesCorrespondingLevels(:,bestMatchesFinerScaleIndicator);
    marginsVector = (-desiredMarginsFromPatchCenter(end):desiredMarginsFromPatchCenter(end));
    marginsTwoDmask = repmat(bsxfun(@plus,abs(marginsVector.'),abs(marginsVector)),[1,1,sum(curScaleIndices2keep),2]);
    marginsTwoDmask(bsxfun(@gt,marginsTwoDmask,shiftdim(desiredMarginsFromPatchCenter(curScaleCorrespondingScales(:,curScaleIndices2keep)).',-2))) = nan;
    marginsTwoDmask(~isnan(marginsTwoDmask)) = 1;
    firstPatchRow = bsxfun(@plus,reshape(firstPatchRow(curScaleIndices2keep),1,1,length(firstPatchRow(curScaleIndices2keep))),marginsVector.');
    firstPatchRow(firstPatchRow<1 | firstPatchRow>largestQueriesMapDimensions(1)) = nan;
    firstPatchCol = bsxfun(@plus,reshape(firstPatchCol(curScaleIndices2keep),1,1,length(firstPatchCol(curScaleIndices2keep))),marginsVector);
    firstPatchCol(firstPatchCol<1 | firstPatchCol>largestQueriesMapDimensions(2)) = nan;
    secondPatchRow = bsxfun(@plus,reshape(secondPatchRow(curScaleIndices2keep),1,1,length(secondPatchRow(curScaleIndices2keep))),marginsVector.');
    secondPatchRow(secondPatchRow<1 | secondPatchRow>largestQueriesMapDimensions(1)) = nan;
    secondPatchCol = bsxfun(@plus,reshape(secondPatchCol(curScaleIndices2keep),1,1,length(secondPatchCol(curScaleIndices2keep))),marginsVector);
    secondPatchCol(secondPatchCol<1 | secondPatchCol>largestQueriesMapDimensions(2)) = nan;
    blockedIndices_firstPatch = sub2ind(largestQueriesMapDimensions,repmat(firstPatchRow,1,length(marginsVector)),repmat(firstPatchCol,length(marginsVector),1));
    blockedIndices_secondPatch = sub2ind(largestQueriesMapDimensions,repmat(secondPatchRow,1,length(marginsVector)),repmat(secondPatchCol,length(marginsVector),1));
    %     Adding inter-pair's patches shifts - adding shifted versions of first
    %     patch's indices, to account for shifts of the first patch relative to the second:
    %     The second patch allways originates from the finer (or equal) scale (lower
    %     scale num). That's why I add the inter-pair shifts to it.
    blockedIndices_secondPatch_withShifts = nan([partialSize(blockedIndices_secondPatch,1:2),...
        maxBlockedInterPatchesPixelsShift(scaleNum)^2,size(blockedIndices_secondPatch,3)]);
    shiftsVect = -floor(maxBlockedInterPatchesPixelsShift(scaleNum)/2):floor(maxBlockedInterPatchesPixelsShift(scaleNum)/2);
    for shiftNum_row = 1:maxBlockedInterPatchesPixelsShift(scaleNum)
        for shiftNum_col = 1:maxBlockedInterPatchesPixelsShift(scaleNum)
            withShifts_curIndices = [max(1,1+shiftsVect(shiftNum_row)),...
                min(size(blockedIndices_secondPatch_withShifts,1),size(blockedIndices_secondPatch_withShifts,1)+shiftsVect(shiftNum_row)),...
                max(1,1+shiftsVect(shiftNum_col)),min(size(blockedIndices_secondPatch_withShifts,1),size(blockedIndices_secondPatch_withShifts,1)+shiftsVect(shiftNum_col))];
            %             The third dimension of blockedIndices_firstPatch_withShifts is going to hold shifted versions of blockedIndices_firstPatch.
            %             Each index in this dimension corresponds to a certain rows-cols shift.
            blockedIndices_secondPatch_withShifts(withShifts_curIndices(1):withShifts_curIndices(2),withShifts_curIndices(3):withShifts_curIndices(4),...
                sub2ind(maxBlockedInterPatchesPixelsShift(scaleNum)*ones(1,2),shiftNum_row,shiftNum_col),:) =...
                blockedIndices_secondPatch((withShifts_curIndices(1):withShifts_curIndices(2))-shiftsVect(shiftNum_row),...
                (withShifts_curIndices(3):withShifts_curIndices(4))-shiftsVect(shiftNum_col),:);
        end
    end
    blockedIndices_firstPatch_allScales = reshape(bsxfun(@times,blockedIndices_firstPatch,...
        marginsTwoDmask(:,:,:,1)),[length(marginsVector).^2,...
        1,sum(curScaleIndices2keep)]);
    blockedIndices_secondPatch_allScales = reshape(bsxfun(@times,blockedIndices_secondPatch_withShifts,...
        reshape(marginsTwoDmask(:,:,:,2),[partialSize(marginsTwoDmask,1:2),1,size(marginsTwoDmask,3)])),[length(marginsVector).^2,...
        size(blockedIndices_secondPatch_withShifts,3),sum(curScaleIndices2keep)]);
    blockedIndices{scaleNum} = sub2ind(prod(largestQueriesMapDimensions)*ones(1,2),repmat(blockedIndices_firstPatch_allScales,...
        [1,size(blockedIndices_secondPatch_allScales,2),1]),blockedIndices_secondPatch_allScales);
    blockedIndices{scaleNum} = reshape(blockedIndices{scaleNum},[prod(partialSize(blockedIndices{scaleNum},1:2)),size(blockedIndices{scaleNum},3)]);
    totalNumOfNonNanBlockedIndices = totalNumOfNonNanBlockedIndices+sum(sum(~isnan(blockedIndices{scaleNum})));
end
clear blockedIndices_firstPatch_withShifts blockedIndices_firstPatch blockedIndices_secondPatch
if size(queriesMapsDimensions,1)>1
    columnedGlobalBlockedIndices = nan(totalNumOfNonNanBlockedIndices,1);
    columnedGlobalBlocked_J = nan(size(columnedGlobalBlockedIndices));  columnedGlobalBlocked_I = nan(size(columnedGlobalBlockedIndices));
    nextFirstIndex = 1;
    for pairNum = 1:length(pairsPriorityOrder)
        curBlockedIndices = blockedIndices{pairCellNumInblockedIndices(pairsPriorityOrder(pairNum))}...
            (:,pairNumMapping2withinScaleValidPairNum(pairsPriorityOrder(pairNum)));
        curBlockedIndices(isnan(curBlockedIndices)) = [];
        columnedGlobalBlockedIndices(nextFirstIndex:nextFirstIndex+length(curBlockedIndices)-1) = curBlockedIndices;
        columnedGlobalBlocked_J(nextFirstIndex:nextFirstIndex+length(curBlockedIndices)-1) = pairNum;
        columnedGlobalBlocked_I(nextFirstIndex:nextFirstIndex+length(curBlockedIndices)-1) = 1:length(curBlockedIndices);
        nextFirstIndex = nextFirstIndex+length(curBlockedIndices);
    end
    clear blockedIndices
    if nextFirstIndex-1~=length(columnedGlobalBlockedIndices);  error('Sanity check failed');   end
else
    globalBlockedIndices = zeros(maxNumOfIndices,length(pairsPriorityOrder));
    globalIndexOfFirstPairInCurScale = 1;
    for scaleNum = 1:size(queriesMapsDimensions,1)
        globalBlockedIndices(1:size(blockedIndices{scaleNum},1),globalIndexOfFirstPairInCurScale+(0:size(blockedIndices{scaleNum},2)-1)) = blockedIndices{scaleNum};
        globalIndexOfFirstPairInCurScale = globalIndexOfFirstPairInCurScale+size(blockedIndices{scaleNum},2);
    end
    nansIncluded2nansExcludedIndicesMapping = nan(size(pairsPriorityOrder));   nansIncluded2nansExcludedIndicesMapping(indices2keep) = 1:sum(indices2keep);
    [columnedGlobalBlocked_I,columnedGlobalBlocked_J,columnedGlobalBlockedIndices] =...
        find(globalBlockedIndices(:,nansIncluded2nansExcludedIndicesMapping(pairsPriorityOrder)));
    columnedGlobalBlocked_I(isnan(columnedGlobalBlockedIndices)) = [];    columnedGlobalBlocked_J(isnan(columnedGlobalBlockedIndices)) = [];
    columnedGlobalBlockedIndices(isnan(columnedGlobalBlockedIndices)) = [];
end
pseudoBlockedIndicesSize = [max(columnedGlobalBlocked_I),max(columnedGlobalBlocked_J)];
validIndicesInglobalBIndices = sub2ind(pseudoBlockedIndicesSize,columnedGlobalBlocked_I,columnedGlobalBlocked_J);
[~,indicesFirstOccurence] = unique(columnedGlobalBlockedIndices);
indicesNotAppearingFirst = setdiff(validIndicesInglobalBIndices,validIndicesInglobalBIndices(indicesFirstOccurence));
[~,indicesNotAppearingFirst_J] = ind2sub(pseudoBlockedIndicesSize,indicesNotAppearingFirst);
invalidIndices = pairsPriorityOrder(unique(indicesNotAppearingFirst_J));
indices2keep(invalidIndices) = false;
if sum(indices2keep)<numOfActuallyDesiredPairs && length(pairsPriorityOrder)<initialPairsPriorityVectorLength
    error('Pairs'' vector was clipped too much...');
end