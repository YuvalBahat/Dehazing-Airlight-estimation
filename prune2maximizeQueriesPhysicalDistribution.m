function indices2keep = prune2maximizeQueriesPhysicalDistribution(queriesIndeces,queryQualityValues,numOfDesiredQueriesOrAllowedQueryOverlap,queriesMapDimensions,samplingWinSize)

if size(queriesIndeces,1)==2; %Dealing with pairs of patches
%     error('No longer supported. Use prune2preventRepeatedPairs instead.');
    indices2keep = true(1,size(queriesIndeces,2));
    indices2keep(any(isnan(queryQualityValues),1)) = false; %Pairs that were a-priori eliminated were indicated by passing a nan as their patchesValues
    if numOfDesiredQueriesOrAllowedQueryOverlap==1 || all(~indices2keep)
        return;
    end
    patchSize = (samplingWinSize+1)/2;
    maxOverlappingPixels = round(numOfDesiredQueriesOrAllowedQueryOverlap*patchSize^2);
    if size(queryQualityValues,1)>1
        ordersByDifferentMeasures = nan(size(queryQualityValues));
        for qualityMeasureNum = 1:size(queryQualityValues,1)
            values2sort = queryQualityValues(qualityMeasureNum,:);  values2sort(isnan(values2sort)) = -realmax;
            [~,ordersByDifferentMeasures(qualityMeasureNum,:)] = sort(values2sort,'descend');
        end
        descendingPairsQualityOrder = mergeMultipleSorts(ordersByDifferentMeasures);
        descendingPairsQualityOrder(any(isnan(queryQualityValues(:,descendingPairsQualityOrder)),1)) = [];
    else
        [~,descendingPairsQualityOrder] = sort(queryQualityValues,'descend');  firstNonNanInd = find(~isnan(queryQualityValues(descendingPairsQualityOrder)),1,'first');
        descendingPairsQualityOrder = descendingPairsQualityOrder(firstNonNanInd:end);
    end
    [firstPatchRow,firstPatchCol] = ind2sub(queriesMapDimensions,queriesIndeces(1,:));
    [secondPatchRow,secondPatchCol] = ind2sub(queriesMapDimensions,queriesIndeces(2,:));
    if 0 %old method, using NN search. Discarding pairs in which one patch shifts in regard to other, because it is considered overlap here.
        pairsSubIndices = [firstPatchRow;firstPatchCol;secondPatchRow;secondPatchCol];
        numOfNeighbors2search = min([80*patchSize,2e3,sum(all(~isnan(queryQualityValues),1))]);
        NNsMat = nan(numOfNeighbors2search,size(queryQualityValues,2));
        nansExcluded2nansIncludedIndecesMapping = find(all(~isnan(queryQualityValues),1));
        [NNsMat(:,all(~isnan(queryQualityValues),1)),~,searchDuration] = ...
            complex_KDSearch(single(pairsSubIndices(:,all(~isnan(queryQualityValues),1))),...%Pairs that were a-priori eliminated were indicated by passing a nan as their patchesValues
            single(pairsSubIndices(:,all(~isnan(queryQualityValues),1))),numOfNeighbors2search,1);
        tocAndPrintWhenLarge(2,'Spatially pruning patches pairs for A estimation',searchDuration);
        NNsMat(:,all(~isnan(queryQualityValues),1)) = nansExcluded2nansIncludedIndecesMapping(NNsMat(:,all(~isnan(queryQualityValues),1)));
        NNsMat(1,:) = [];
        for pairNum = descendingPairsQualityOrder
            if indices2keep(pairNum)
                overlapInAxis = reshape(patchSize-abs(bsxfun(@minus,pairsSubIndices(:,pairNum),pairsSubIndices(:,NNsMat(:,pairNum)))),2,2,[]);
                overlapInAxis(overlapInAxis<0) = 0;
                tooCloseNeighbors = all(prod(overlapInAxis)/maxOverlappingPixels>1,2);
                indices2keep(NNsMat(tooCloseNeighbors,pairNum)) = false;
            end
        end
    else %New method, using indeces in the numel(image)x numel(image) domain, and unique to find repetitions. Allowing shifts in one patch of the pair
        desiredMarginsFromPatchCenter = ceil((1-numOfDesiredQueriesOrAllowedQueryOverlap)*patchSize-1);
        marginsVector = (-desiredMarginsFromPatchCenter:desiredMarginsFromPatchCenter);
        firstPatchRow = bsxfun(@plus,shiftdim(firstPatchRow(indices2keep),-1),marginsVector.');  firstPatchRow(firstPatchRow<1 | firstPatchRow>queriesMapDimensions(1)) = nan;
        firstPatchCol = bsxfun(@plus,shiftdim(firstPatchCol(indices2keep),-1),marginsVector);  firstPatchCol(firstPatchCol<1 | firstPatchCol>queriesMapDimensions(2)) = nan;
        secondPatchRow = bsxfun(@plus,shiftdim(secondPatchRow(indices2keep),-1),marginsVector.');  secondPatchRow(secondPatchRow<1 | secondPatchRow>queriesMapDimensions(1)) = nan;
        secondPatchCol = bsxfun(@plus,shiftdim(secondPatchCol(indices2keep),-1),marginsVector);  secondPatchCol(secondPatchCol<1 | secondPatchCol>queriesMapDimensions(2)) = nan;
        blockedIndices_firstPatch = reshape(sub2ind(queriesMapDimensions,repmat(firstPatchRow,1,length(marginsVector)),repmat(firstPatchCol,length(marginsVector),1)),...
            [length(marginsVector).^2,size(descendingPairsQualityOrder,2)]);
        blockedIndices_secondPatch = reshape(sub2ind(queriesMapDimensions,repmat(secondPatchRow,1,length(marginsVector)),repmat(secondPatchCol,length(marginsVector),1)),...
            [length(marginsVector).^2,size(descendingPairsQualityOrder,2)]);
        blockedIndeces = sub2ind(prod(queriesMapDimensions)*ones(1,2),blockedIndices_secondPatch,blockedIndices_firstPatch);
        nansIncluded2nansExcludedIndecesMapping = nan(size(descendingPairsQualityOrder));   nansIncluded2nansExcludedIndecesMapping(indices2keep) = 1:sum(indices2keep);
        [~,indecesFirstOccurence] = unique(reshape(blockedIndeces(:,nansIncluded2nansExcludedIndecesMapping(descendingPairsQualityOrder)),[],1));
        indecesFirstOccurenceMat = false(size(blockedIndeces)); indecesFirstOccurenceMat(indecesFirstOccurence) = true;
        indices2keep(descendingPairsQualityOrder) = all(indecesFirstOccurenceMat,1);
%         patchPairsMat = spalloc(prod(queriesMapDimensions),prod(queriesMapDimensions),numel(blockedIndeces));
%         %         blockedMatrixIndeces = nan(1,numel(blockedIndeces));    blockedMatrixIndecesVectIndex = 1;
%         pairsCounter = 0;
%         for pairNum = descendingPairsQualityOrder
%             if indices2keep(pairNum)
%                 relevantPairIndeces = blockedIndeces(:,nansIncluded2nansExcludedIndecesMapping(pairNum)); relevantPairIndeces(isnan(relevantPairIndeces)) = [];
%                 if any(patchPairsMat(relevantPairIndeces))
%                     indices2keep(pairNum) = false;
%                 else
%                     pairsCounter = pairsCounter+1;
%                     tic
%                     patchPairsMat(relevantPairIndeces) = true;
%                     time2assign(pairsCounter) = toc;
%                 end
%             end
%         end
    end
else
    if size(queryQualityValues,1)>1;    error('Unsupported');   end;
    desiredQueriesPortion = numOfDesiredQueriesOrAllowedQueryOverlap/prod(queriesMapDimensions);
    if desiredQueriesPortion>=1
        indices2keep = queriesIndeces;
        return;
    end
    
    queriesLocationsMap = zeros(queriesMapDimensions);    queriesValuesMap = nan(size(queriesLocationsMap));    %valThreshold = min(patchesValues);
    queriesLocationsMap(queriesIndeces) = 1;  queriesValuesMap(queriesIndeces) = queryQualityValues;
    queriesDensity = conv2(padarray(queriesLocationsMap,floor(samplingWinSize/2)*ones(1,2),'replicate','both'),1/samplingWinSize^2*ones(samplingWinSize),'valid');
    
    sortedDensity = sort(queriesDensity(:));
    if sortedDensity(1)>desiredQueriesPortion
        densityThreshold = desiredQueriesPortion;
    else
        leftQueries = prod(queriesMapDimensions):-1:1;
        densityThreshold = sortedDensity(find((cumsum(sortedDensity)+leftQueries.'.*sortedDensity)/prod(queriesMapDimensions)>desiredQueriesPortion,1,'first'));
    end
    smoothedPatchesValuesMap = smooth2a(queriesValuesMap,floor(samplingWinSize/2),floor(samplingWinSize/2));
    threshold4stochasticPruning = min(densityThreshold./queriesDensity.*queriesValuesMap./smoothedPatchesValuesMap/...
        (mean(mean(queriesValuesMap(smoothedPatchesValuesMap>0)./smoothedPatchesValuesMap(smoothedPatchesValuesMap>0)))),1);
    queriesLocationsMap = (rand(queriesMapDimensions)<=threshold4stochasticPruning).*queriesLocationsMap;
    indices2keep = logical(queriesLocationsMap(queriesIndeces)).';
end
end

function extendedIndices = returnExtendedPatchesIndices(orgPatchIndices,patchesMapDimensions,patchSize)
halfPatchSize = floor(patchSize/2);
[firstPatchRow,firstPatchCol] = ind2sub(patchesMapDimensions,orgPatchIndices);
firstPatchRow = bsxfun(@plus,firstPatchRow,(-halfPatchSize:1:halfPatchSize).');
firstPatchCol = bsxfun(@plus,firstPatchCol,(-halfPatchSize:1:halfPatchSize).');
extendedPatchIndicesCols = nan(size(firstPatchRow,1)^2,size(orgPatchIndices,2));
extendedPatchIndicesRows = nan(size(firstPatchRow,1)^2,size(orgPatchIndices,2));
for i = 1:size(extendedPatchIndicesRows,2)
    [curCols,curRows] = meshgrid(firstPatchCol(:,i),firstPatchRow(:,i));
    extendedPatchIndicesCols(:,i) = curCols(:);
    extendedPatchIndicesRows(:,i) = curRows(:);
end
extendedPatchIndicesRows(extendedPatchIndicesRows>patchesMapDimensions(1) | extendedPatchIndicesRows<1) = nan;
extendedPatchIndicesCols(extendedPatchIndicesCols>patchesMapDimensions(2) | extendedPatchIndicesCols<1) = nan;
extendedIndices = sub2ind(patchesMapDimensions,extendedPatchIndicesRows,extendedPatchIndicesCols);
end