function pairsValidity = prioritizeNoOverlapAndPrunePairs(values,originalFirstPatches,bestMatches_perScalePatchIndeces,bestMatchesCorrespondingScales,...
    rescaledPatchesMapDimensions,overlapPruningPrioritizingMethod,maxAllowedPairsOverlap,patchSize,desiredNumOfPairsRange,minNormalizedCorrelation,...
    pruneOverlap)

if ~pruneOverlap && (sum(~isnan(values.overallNNsValidityVals))<=desiredNumOfPairsRange(2))
    pairsPriorityOrder = find(~isnan(values.overallNNsValidityVals));
else
    switch overlapPruningPrioritizingMethod
        case 'reoccurrenceScore'
            values.reoccurrenceVals(isnan(values.overallNNsValidityVals)) = nan;
            [~,pairsPriorityOrder] = sort(values.reoccurrenceVals);
            pairsPriorityOrder(isnan(values.reoccurrenceVals(pairsPriorityOrder))) = [];
        case 'complexity_meansVariety_tRatios'
            firstPatchMeans = reshape(mean(originalFirstPatches),partialSize(originalFirstPatches,3:4));
            firstPatchMeans(isnan(values.overallNNsValidityVals),:) = nan;
            prioritizingVals.correspondingShiftedImageNC = max(values.bestMatchesCorrespondingShiftsNC);
            prioritizingVals.tRatios = values.tRatios;
            pairsPriorityOrder = returnOrderMaximizingPairsMeansVariaty(firstPatchMeans,values.overallNNsValidityVals,100,prioritizingVals);
    end
end
if length(desiredNumOfPairsRange)<2;    desiredNumOfPairsRange(2) = size(bestMatches_perScalePatchIndeces,2);   end;
if isempty(pairsPriorityOrder); error('All pairs discarded');   end;
if pruneOverlap
    pairsValidity = prune2preventRepeatedPairs(bestMatches_perScalePatchIndeces,bestMatchesCorrespondingScales,pairsPriorityOrder,...
        maxAllowedPairsOverlap,rescaledPatchesMapDimensions,patchSize,0,desiredNumOfPairsRange(2));
else
    pairsValidity = ~isnan(values.overallNNsValidityVals);
end
pairsPriorityOrder(~pairsValidity(pairsPriorityOrder)) = [];
if any(pairsValidity)
    numOfPairs2use = sum(values.overallNNsValidityVals(pairsPriorityOrder)>=minNormalizedCorrelation);
    if numOfPairs2use>=desiredNumOfPairsRange(1)
        pairsPriorityOrder = pairsPriorityOrder(values.overallNNsValidityVals(pairsPriorityOrder)>=minNormalizedCorrelation);
        numOfPairs2use = min(length(pairsPriorityOrder),desiredNumOfPairsRange(2));
    else
        NCthreshold = sort(values.overallNNsValidityVals(pairsPriorityOrder),'descend');
        NCthreshold = NCthreshold(min(length(NCthreshold),desiredNumOfPairsRange(1)));
        pairsPriorityOrder = pairsPriorityOrder(values.overallNNsValidityVals(pairsPriorityOrder)>=NCthreshold);
        numOfPairs2use = length(pairsPriorityOrder);
    end
    pairsValidity = false(size(pairsValidity));   pairsValidity(pairsPriorityOrder(1:numOfPairs2use)) = true;
end