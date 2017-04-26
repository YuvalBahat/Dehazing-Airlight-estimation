function [t2DividedByt1,orderedFactorizationErrors] = computePatchesRatio(NNpatches,tRatioCalculation,patchesGroupNums,patchesNC,NCthreshold)
separateT4eachChannel = false;
if separateT4eachChannel
    switch tRatioCalculation
        case 'norm2Ratios'
            t2DividedByt1 = sqrt(bsxfun(@rdivide,sum((NNpatches(:,2:end,:,:)).^2,1),sum((NNpatches(:,1,:,:)).^2,1)));
        case 'normalizedInnerProduct'
            t2DividedByt1 = sum(bsxfun(@times,NNpatches(:,1,:,:),NNpatches(:,2:end,:,:)))./sum((NNpatches(:,1,:,:)).^2);
        case 'norm1Ratios'
            t2DividedByt1 = bsxfun(@rdivide,sum(abs(NNpatches(:,2:end,:,:))),sum(abs(NNpatches(:,1,:,:))));
    end
else
    switch tRatioCalculation
        case 'norm1Ratios'
            t2DividedByt1 = sum(sum(abs(NNpatches(:,2:end,:,:)),1),4)./sum(sum(abs(NNpatches(:,1,:,:)),1),4);
        case 'norm2Ratios'
            t2DividedByt1 = sqrt(sum(sum((NNpatches(:,2:end,:,:)).^2,1),4)./sum(sum((NNpatches(:,1,:,:)).^2,1),4));
        case 'factorization'
            [indexInMatrix,orderOfpatches2Use,maxGroupRecurrence,numOfUniqueGroups] = returnGatherByCommonIdentifiersStructer(patchesGroupNums,[],true);
            patchesNC(patchesNC<0.01) = 0.01;
            patchesAssumedNoiseVar = nan(maxGroupRecurrence,numOfUniqueGroups);   patchesAssumedNoiseVar(indexInMatrix) = patchesNC(orderOfpatches2Use);
            NCwiseValidity = patchesAssumedNoiseVar;    NCwiseValidity = cat(1,true(1,numOfUniqueGroups),NCwiseValidity>=NCthreshold);
            patchesAssumedNoiseVar = shiftdim(1./patchesAssumedNoiseVar,-1);
            groupedNNpatches = nan(prod(partialSize(NNpatches,[1,4])),maxGroupRecurrence+1,numOfUniqueGroups);
            [indexInMatrix_row,indexInMatrix_col] = ind2sub([maxGroupRecurrence,numOfUniqueGroups],indexInMatrix);
            indexInMatrix = sub2ind(size(groupedNNpatches),repmat((1:prod(partialSize(NNpatches,[1,4]))).',[1,length(indexInMatrix_row)]),...
                repmat(indexInMatrix_row.'+1,[prod(partialSize(NNpatches,[1,4])),1]),repmat(indexInMatrix_col.',[prod(partialSize(NNpatches,[1,4])),1]));
            NNpatches = reshape(permute(NNpatches,[1,4,2,3]),[prod(partialSize(NNpatches,[1,4])),2,size(NNpatches,3)]);
            [~,indicesOfGroupsFirstPatches] = unique(patchesGroupNums,'stable');
            groupedNNpatches(:,1,:) = NNpatches(:,1,indicesOfGroupsFirstPatches);
            groupedNNpatches(indexInMatrix) = NNpatches(:,2,orderOfpatches2Use);
            [~,tDividedByAlpha,factorizationErr] = alternativeRank1MatrixFactorization(groupedNNpatches,cat(2,ones(1,1,numOfUniqueGroups),patchesAssumedNoiseVar),...
                shiftdim(all(~isnan(groupedNNpatches),1),1) & NCwiseValidity);
            tDividedByAlpha = shiftdim(tDividedByAlpha,1);
            t2DividedByt1_preOrder = bsxfun(@rdivide,tDividedByAlpha(2:end,:),tDividedByAlpha(1,:));
            if nargout>1 %for debugging
                %%
                processedFactorizationErr = shiftdim(factorizationErr(1,2:end,:),1); processedFactorizationErr = processedFactorizationErr(~isnan(t2DividedByt1_preOrder));
                orderedFactorizationErrors = nan(size(patchesNC));  
                orderedFactorizationErrors(orderOfpatches2Use(patchesNC(orderOfpatches2Use)>=NCthreshold)) = processedFactorizationErr;
                if false
                    %%
                    figure;
                    scatter(patchesNC(patchesNC>=NCthreshold),orderedFactorizationErrors(patchesNC>=NCthreshold));   maximizeFigure;
                    figure;
                    multiD21Dhist(processedFactorizationErr,0);
                end
            end
            t2DividedByt1(orderOfpatches2Use(patchesNC(orderOfpatches2Use)>=NCthreshold)) = t2DividedByt1_preOrder(~isnan(t2DividedByt1_preOrder));
            t2DividedByt1(patchesNC<NCthreshold) = nan;
            t2DividedByt1 = reshape(t2DividedByt1,[1,1,length(patchesNC)]);
        otherwise
            error('Not implemented yet');
    end
end