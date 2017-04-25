function [estimatedAirlight,individualPairAirlight,individualPairAvar] = estimateAirlightColor(originalChosenPatches,DClessChosenPatches,...
    estimatedNoiseVarOrAguess,method,t1Estimates)

Uvector = DClessChosenPatches(:,1,:,:)-DClessChosenPatches(:,2,:,:);
Vvector = originalChosenPatches(:,2,:,:).*DClessChosenPatches(:,1,:,:)-...
    originalChosenPatches(:,1,:,:).*DClessChosenPatches(:,2,:,:);
switch method
    case 'directWithoutNoise'
        estimatedNoiseVarOrAguess = estimatedNoiseVarOrAguess{1};
        Uvector = bsxfun(@times,reshape(1./estimatedNoiseVarOrAguess,1,1,[]),Uvector);
        Vvector = bsxfun(@times,reshape(1./estimatedNoiseVarOrAguess,1,1,[]),Vvector);
        curAirlightVectEst_nominator = reshape(nanmean(Uvector.*Vvector,1),partialSize(Uvector,3:4));
        curAirlightVectEst_denominator = reshape(nanmean(Uvector.^2,1),partialSize(Uvector,3:4));
        individualPairAirlight = curAirlightVectEst_nominator./curAirlightVectEst_denominator;
        estimatedAirlight = sum(curAirlightVectEst_nominator,1)./sum(curAirlightVectEst_denominator,1);
        if nargout>2
            DClessPatchNorms = sqrt(sum(DClessChosenPatches.^2,1));
            individualPairAvar = max(shiftdim(var(bsxfun(@rdivide,bsxfun(@times,DClessPatchNorms(1,2,:,:),...
                originalChosenPatches(:,1,:,:))-bsxfun(@times,DClessPatchNorms(1,1,:,:),...
                originalChosenPatches(:,2,:,:)),DClessPatchNorms(1,2,:,:)-DClessPatchNorms(1,1,:,:)),0,1),2),[],2);
        end
    case 'directWithNoise'
        maxIterations = 100;
        I2mean = shiftdim(nanmean(originalChosenPatches(:,2,:,:),1),2);
        previousAirlightEst = zeros(1,3);
        tempAirlight = [];  tempInvalidity = [];    tempAminusI2Mean = [];
        for i=1:maxIterations
            if i>maxIterations-1
                error('Too many airlight noise model estimation iterations');
            end
            % Creating the covariance matrices:
            if i==1
                noiseVariances = 0.01*ones(size(I2mean,1),1);
            else
                if i==2
                    currentNoiseVariances = [];
                end
                if isempty(t1Estimates)
                    tLB = reshape(nanmax(returnTransmissionLbasedLowerBound(permute(originalChosenPatches,[3,4,1,2]),estimatedAirlight),[],3),...
                        partialSize(originalChosenPatches,[3,2]));
                else
                    tLB = t1Estimates;
                end
                noiseVariances = -diff(tLB,1,2).*(tLB(:,1)./max(0.01,tLB(:,2))-1);
                currentNoiseVariances = cat(3,currentNoiseVariances,noiseVariances);
                noiseVariances = (1./mean(currentNoiseVariances(:,:,(max(1,floor(i/2)):i-1)),3)).^2;
            end
            individualPairAirlight = shiftdim(nansum(Uvector.*Vvector,1),2)./shiftdim(nansum(Uvector.^2,1),2);
            weights = 1./noiseVariances;
            estimatedAirlight = sum(bsxfun(@times,weights,individualPairAirlight),1)./sum(weights,1);
            tempAirlight = [tempAirlight;estimatedAirlight];
            if norm(previousAirlightEst-estimatedAirlight)<1e-3 %|| i>maxIterations
                %                     i
                break
            end
            previousAirlightEst = estimatedAirlight;
        end
end