function rescaledLocations = rescalePatchMapLocations2HighestLevel(originalLocations,correspondingMapLevels,queriesMapsDimensions,patchSizeUsed4PatchMaps,...
    scaleNum2choose)
rescaledLocations = nan(size(originalLocations));
if exist('scaleNum2choose','var')
    chosenDimsIndex = scaleNum2choose;
else
    [~,chosenDimsIndex] = max(queriesMapsDimensions(:,1));
end
largestQueriesMapDimensions = queriesMapsDimensions(chosenDimsIndex,:)+patchSizeUsed4PatchMaps-1;
downscalingFactors = largestQueriesMapDimensions(1)./(queriesMapsDimensions(:,1)+patchSizeUsed4PatchMaps-1);
for scaleNum = 1:size(queriesMapsDimensions,1)
    curScaleLocationsIndicator = correspondingMapLevels==scaleNum;
    if any(originalLocations(curScaleLocationsIndicator(:))>prod(queriesMapsDimensions(scaleNum,:)));   error('Patch map index exceeds map size');  end;
    [rows,cols] = ind2sub(queriesMapsDimensions(scaleNum,:),originalLocations(curScaleLocationsIndicator));
    if downscalingFactors(scaleNum)~=1 % Transforming row and column's patch indices to highest level size image row & col patch indices:
        rows = round((rows+floor(patchSizeUsed4PatchMaps/2))*downscalingFactors(scaleNum)-floor(patchSizeUsed4PatchMaps/2));
        cols = round((cols+floor(patchSizeUsed4PatchMaps/2))*downscalingFactors(scaleNum)-floor(patchSizeUsed4PatchMaps/2));
%         if downscalingFactors(scaleNum)<1
            rows = min(queriesMapsDimensions(chosenDimsIndex,1),max(1,rows));
            cols = min(queriesMapsDimensions(chosenDimsIndex,2),max(1,cols));
%         end
    end
    rescaledLocations(curScaleLocationsIndicator) = sub2ind(queriesMapsDimensions(chosenDimsIndex,:),rows,cols);
end