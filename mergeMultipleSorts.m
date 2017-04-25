function indecesMergedOrder = mergeMultipleSorts(sortedIndeesOrder)
% Merging differently sorted indeces of an array, yeilding a new indeces
% order. The first index in the new order is the one that is first to
% appear in ALL differently ordered input orders. The last one is the that
% is the last to appear in all different sorts.
% If 2 (or more) indices switch places in the different rankings, they will
% be sorted according to the ranking appearing in a greater row. For
% example:
% mergeMultipleSorts([1,2;2,1]) = [2,1], while
% mergeMultipleSorts([2,1;1,2]) = [1,2]. This means the more improtant
% ranking should be put to a greater row index in sortedIndeesOrder.
% Inputs:
% sortedIndeesOrder - an MxNxL array: M different orders of indeces 1:N. L
% different (independent) orderings
% Outputs:
% indecesMergedOrder - The resulting order of indeces.
C = nan(partialSize(sortedIndeesOrder,2:3));    IA = C;
for i=1:size(sortedIndeesOrder,3)
    [C(:,i),IA(:,i),~] = unique(reshape(sortedIndeesOrder(:,:,i),[],1),'last');
end
if ispc
    [~,memStatus] = memory; IA_Properties = whos('IA');
end
if ~ispc || memStatus.PhysicalMemory.Available<(IA_Properties.bytes)
    numOfIterations = ceil(size(IA,2)/1e3);
    indeces2takeFromC = nan(size(IA));
else
    numOfIterations = 1;
end
for iter = 1:numOfIterations
    if numOfIterations>1;   thisIterQueryIndices = (iter-1)*1e3+1:min(iter*1e3,size(IA,2));   else   thisIterQueryIndices = 1:size(IA,2);   end
    [~,indeces2takeFromC(:,thisIterQueryIndices)] = sort(IA(:,thisIterQueryIndices),'ascend');
end
indecesMergedOrder = C(bsxfun(@plus,indeces2takeFromC,0:size(C,1):numel(C)-1)).';