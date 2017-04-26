function resizedImage = imresize2fitDesiredSize(noHazeIm,desiredImSize,sizeCalculation,resizeMethod)
switch sizeCalculation
    case 'minDim'
        [currentSize,relevantDim] = min(partialSize(noHazeIm,1:2));
        newImSize(relevantDim) = desiredImSize;
        newImSize(one2Nmod(relevantDim-1,2)) = round(desiredImSize/currentSize*size(noHazeIm,one2Nmod(relevantDim-1,2)));
%         relativeDownscaleFactor = ceil(min(orgImSize)/scalesLegend(downscaleFactor));
    case 'maxDim'
        [currentSize,relevantDim] = max(partialSize(noHazeIm,1:2));
        newImSize(relevantDim) = desiredImSize;
        newImSize(one2Nmod(relevantDim-1,2)) = round(desiredImSize/currentSize*size(noHazeIm,one2Nmod(relevantDim-1,2)));
        error('code not checked yet');
    case 'area'
        currentArea = prod(partialSize(noHazeIm,1:2));
        newImSize = round(sqrt(desiredImSize/currentArea)*partialSize(noHazeIm,1:2));
    case 'averageSize'
        currentSize = geomean(partialSize(noHazeIm,1:2));
        newImSize = round(sqrt(desiredImSize/currentSize)*partialSize(noHazeIm,1:2));
    otherwise
        error('Unsupported');
end
if ~exist('resizeMethod','var')
    resizeMethod = 'bicubic';
end
resizedImage = imresize(noHazeIm,newImSize,resizeMethod);