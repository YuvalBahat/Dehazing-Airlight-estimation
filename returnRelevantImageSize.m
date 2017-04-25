function relevantImSize = returnRelevantImageSize(orgImSize,method)
switch method
    case 'minDim'
        relevantImSize = min(orgImSize);
%         relativeDownscaleFactor = ceil(min(orgImSize)/scalesLegend(downscaleFactor));
    case 'maxDim'
        relevantImSize = max(orgImSize);
    case 'averageSize'
        relevantImSize = geomean(orgImSize);
    case 'area'
        relevantImSize = prod(orgImSize);
    otherwise
        error('Unsupported');
end