function [foldIndeces,patchesAssignment,initialTransmissionMapGuess,foldsAssignmentMap] = assignPatches2Folds(patchIndices,patchesMapDimension,numOfFolds,foldsDivisionMethod,...
    patchesSTDvalues,image,estimatedAirlight,patchSize,foldsAssignmentMap)
foldIndeces = cell(numOfFolds,1);patchesAssignment = nan(size(patchIndices));
switch foldsDivisionMethod
    case 'horizontalStripes'
        [~,initialTransmissionMapGuess] = meshgrid(1:patchesMapDimension(2),1:patchesMapDimension(1));
        initialTransmissionMapGuess = -1*initialTransmissionMapGuess;
        [rowNums,~] = ind2sub(patchesMapDimension,patchIndices);
        foldsRowBoundries = round(linspace(0,patchesMapDimension(1),numOfFolds+1));
        for foldNum = 1:numOfFolds
            foldIndeces{foldNum} = find(rowNums>foldsRowBoundries(foldNum) & rowNums<=foldsRowBoundries(foldNum+1));
            patchesAssignment(foldIndeces{foldNum}) = foldNum;
        end
    case {'brightness','variance','airlightBasedCombination','combinedLowerBound','brightnessVarCombination','tLB'}
        if strcmp(foldsDivisionMethod,'brightnessVarCombination')
            [~,~,V] = rgb2hsv(image);
            brightnessImage = 1-conv2(V,ones(patchSize)/(patchSize^2),'valid');
            varsImage = reshape(patchesSTDvalues,partialSize(image,1:2)-patchSize*ones(1,2)+1);
            initialTransmissionMapGuess = max(cat(3,(brightnessImage-min(brightnessImage(:)))/(max(brightnessImage(:))-min(brightnessImage(:))),...
                (varsImage-min(varsImage(:)))/(max(varsImage(:))-min(varsImage(:)))),[],3);
        elseif strcmp(foldsDivisionMethod,'combinedLowerBound')
%             initialTransmissionMapGuess =...
%                 cat(3,reshape(patchesSTDvalues,patchesMapDimension),...
%                 max(bsxfun(@rdivide,bsxfun(@minus,image,reshape(estimatedAirlight,1,1,3)),1-reshape(estimatedAirlight,1,1,3)),[],3),...
%                 max(bsxfun(@rdivide,bsxfun(@minus,reshape(estimatedAirlight,1,1,3),image),reshape(estimatedAirlight,1,1,3)),[],3));
            initialTransmissionMapGuess = max(...
                cat(4,reshape(patchesSTDvalues,[patchesMapDimension,size(patchesSTDvalues,3)]),...
                bsxfun(@rdivide,bsxfun(@minus,image,reshape(estimatedAirlight,1,1,3)),1-reshape(estimatedAirlight,1,1,3)),...
                bsxfun(@rdivide,bsxfun(@minus,reshape(estimatedAirlight,1,1,3),image),reshape(estimatedAirlight,1,1,3))),[],4);
        elseif strcmp(foldsDivisionMethod,'variance')
            initialTransmissionMapGuess = reshape(patchesSTDvalues,[patchesMapDimension,3]);
%             initialDepthMapGuess = -1*log(initialDepthMapGuess-min(initialDepthMapGuess(:))+.1);
        elseif strcmp(foldsDivisionMethod,'tLB')
            initialTransmissionMapGuess = reshape(patchesSTDvalues,patchesMapDimension);
        elseif strcmp(foldsDivisionMethod,'brightness')% Brightness:
            initialTransmissionMapGuess = -1*rgb2gray(image);
%             initialDepthMapGuess = -1*log(initialDepthMapGuess-min(initialDepthMapGuess(:))+.1);
            initialTransmissionMapGuess = initialTransmissionMapGuess-min(initialTransmissionMapGuess(:));
        else %airlightBasedCombination
            STDimage = reshape(patchesSTDvalues,patchesMapDimension);
            STDimage = ((STDimage)-min((STDimage(:)))+.1);
            patchSize = size(image,1)-patchesMapDimension(1)+1;
            patchMeanImage = conv2_3D(image,1/(patchSize^2)*ones(patchSize),'valid');
            patchMeanImage = rgb2gray(bsxfun(@minus,reshape(estimatedAirlight,1,1,3),patchMeanImage));
            patchMeanImage = (abs(patchMeanImage)-min(abs(patchMeanImage(:)))+.1);
%             initialDepthMapGuess = -1*max(log(patchMeanImage),log(STDimage));
            initialTransmissionMapGuess = max(patchMeanImage,STDimage);
        end
        if nargin<9 %foldsAssignmentMap is not given
            averagingFilter = fspecial('average',100);
            smoothedInitialTransmissionMapGuess = reshape(imfilter(initialTransmissionMapGuess(:,:,1),averagingFilter,'replicate'),1,[]);
            relevantPatchesImage = smoothedInitialTransmissionMapGuess(patchIndices);
            [histVals,binCenters] = multiD21Dhist(relevantPatchesImage,[],'norm');
            [foldsEdges ,~] = find(diff(bsxfun(@minus,histVals,linspace(0,1,numOfFolds+1))>0));
            foldsEdges = [-inf,binCenters(foldsEdges),inf];
            [foldsAssignmentMap,~] = find(bsxfun(@gt,smoothedInitialTransmissionMapGuess,foldsEdges(1:end-1).') & bsxfun(@le,smoothedInitialTransmissionMapGuess,foldsEdges(2:end).'));
            foldsAssignmentMap = reshape(foldsAssignmentMap,size(initialTransmissionMapGuess(:,:,1)));
        end
        patchesAssignment = foldsAssignmentMap(patchIndices);
        %             [patchesAssignment,~] = find(bsxfun(@gt,relevantPatchesImage,foldsEdges(1:end-1).') & bsxfun(@le,relevantPatchesImage,foldsEdges(2:end).'));
        %         previousFoldThreshold = -intmax;
        for foldNum = 1:numOfFolds
            %             curFoldThreshold = binCenters(find(histVals>foldNum/numOfFolds,1,'first'));
            foldIndeces{foldNum} = find(patchesAssignment==foldNum);
            %             previousFoldThreshold = curFoldThreshold;
            %             patchesAssignment(foldIndeces{foldNum}) = foldNum;
        end
        %         foldIndeces{numOfFolds} = find(relevantPatchesImage>curFoldThreshold);
        %         patchesAssignment(foldIndeces{numOfFolds}) = numOfFolds;
    otherwise
        error('Unsupported')
end