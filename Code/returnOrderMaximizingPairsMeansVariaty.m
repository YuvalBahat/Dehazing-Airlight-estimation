function priorityOrder = returnOrderMaximizingPairsMeansVariaty(firstPatchMeans,NCgrades,histGridDenstiy,prioritizingVals)
% Ordering pairs to increase variaty in participating pairs. Dividing the
% 3D space of mean first patch color into histGridDenstiy bins in each
% dimension (total of histGridDenstiy^3 bins). Then ordering pairs by
% taking one pair from each bin, then second pair from each bin, and so on.
% order of bins to take from:
binsOrder = 'NCofFirstPair';%'binSize','NCofFirstPair'

[~,NCorder] = sort(NCgrades,'descend');  NCorder = [NCorder(sum(isnan(NCgrades))+1:end),NCorder(1:sum(isnan(NCgrades)))];
edgesVect = linspace(0,1,histGridDenstiy+1);
[count,~,~,loc] = histcn(firstPatchMeans,edgesVect,edgesVect,edgesVect);
nonEmptyBinsIndicator = count>0;
pairsNonEmptyBinsAssignment = nan(size(loc,1),1);
pairsNonEmptyBinsAssignment(loc(:,1)>0) = sub2ind(size(count),loc(loc(:,1)>0,1),loc(loc(:,1)>0,2),loc(loc(:,1)>0,3));  
switch binsOrder
    case 'NCofFirstPair'
        if nargin>3
            [maxNCperBin,correspondingIndices] = calcMaxNCperBin(pairsNonEmptyBinsAssignment,NCgrades);
            if isfield(prioritizingVals,'correspondingShiftedImageNC')
                [~,nonEmptyBinsOrder] = sort(maxNCperBin.'./prioritizingVals.correspondingShiftedImageNC(correspondingIndices),'descend');
            else
                [~,nonEmptyBinsOrder] = sort(maxNCperBin,'descend');
            end
            if isfield(prioritizingVals,'tRatios')
                [~,tRatiosOrder] = sort(abs(log(prioritizingVals.tRatios(correspondingIndices))),'descend');
                nonEmptyBinsOrder = mergeMultipleSorts([nonEmptyBinsOrder;tRatiosOrder]);
            end
            if isfield(prioritizingVals,'reoccrrenceScores')
                [~,reoccurrenceOrder] = sort(abs(log(prioritizingVals.reoccrrenceScores(correspondingIndices))),'ascend');
                nonEmptyBinsOrder = mergeMultipleSorts([nonEmptyBinsOrder;reoccurrenceOrder]);
            end
        else
            [~,nonEmptyBinsOrder] = sort(calcMaxNCperBin(pairsNonEmptyBinsAssignment,NCgrades),'descend');  
        end
        [~,nonEmptyBinsOrder] = sort(nonEmptyBinsOrder);
    case 'binSize'
        [~,nonEmptyBinsOrder] = sort(count(nonEmptyBinsIndicator(:)),'descend'); [~,nonEmptyBinsOrder] = sort(nonEmptyBinsOrder);
end
binIndices2sortedNonEmptyMapping = nan(size(count));    binIndices2sortedNonEmptyMapping(nonEmptyBinsIndicator) = nonEmptyBinsOrder;
pairsNonEmptyBinsAssignment(loc(:,1)>0) = binIndices2sortedNonEmptyMapping(pairsNonEmptyBinsAssignment(loc(:,1)>0));
[indexInMatrix,orderedNeighborLocations,maxRecurrence] =...
    returnGatherByCommonIdentifiersStructer(pairsNonEmptyBinsAssignment(NCorder),max(pairsNonEmptyBinsAssignment));
[rows,cols] = ind2sub([maxRecurrence,max(pairsNonEmptyBinsAssignment)],indexInMatrix);  
indexInMatrix = sub2ind([max(pairsNonEmptyBinsAssignment),maxRecurrence],cols,rows);
[~,orderOfMatIndices] = sort(indexInMatrix);
priorityOrder = NCorder(orderedNeighborLocations(orderOfMatIndices));
end

function [maxNCperBin,correspondingIndices] = calcMaxNCperBin(pairsNonEmptyBinsAssignment,NCgrades)
noNan2WithNanMapping = find(~isnan(pairsNonEmptyBinsAssignment));
pairsNonEmptyBinsAssignment = pairsNonEmptyBinsAssignment(~isnan(pairsNonEmptyBinsAssignment)); NCgrades = NCgrades(~isnan(NCgrades));
[~,~,uniqueBinNums] = unique(pairsNonEmptyBinsAssignment);
maxNCperBin = zeros(max(uniqueBinNums),1);  correspondingIndices = nan(size(maxNCperBin));
for i=1:length(pairsNonEmptyBinsAssignment)
    if maxNCperBin(uniqueBinNums(i))<NCgrades(i)
        maxNCperBin(uniqueBinNums(i)) = NCgrades(i);
        correspondingIndices(uniqueBinNums(i)) = noNan2WithNanMapping(i);
    end
end
end