noMoreAlgorithmSteps = false;
if ~isempty(parameters.downscaledSize.multiScale)
    if parameters.downscaledSize.multiScale(1)>1
        if algorithmStepNum<=length(parameters.downscaledSize.multiScale)
            downscaleSize = parameters.downscaledSize.multiScale(algorithmStepNum);
        else
            noMoreAlgorithmSteps = true;
        end            
    else
        if downscaleSize==parameters.downscaledSize.multiScale(2) || algorithmStepNum>=parameters.downscaledSize.maxNumOfScales
            noMoreAlgorithmSteps = true;
        else
            downscaleSize = round(parameters.downscaledSize.multiScale(1)*downscaleSize);
            if downscaleSize<parameters.downscaledSize.multiScale(2)/sqrt(parameters.downscaledSize.multiScale(1))
                downscaleSize = parameters.downscaledSize.multiScale(2);
            end
        end
    end
    tempAlgorithmStepNum = algorithmStepNum;
    while downscaleSize>=returnRelevantImageSize(orgImSize,parameters.downscaledSize.approach)
        tempAlgorithmStepNum = tempAlgorithmStepNum+1;
        scalesFirstImPatchIndOffset(tempAlgorithmStepNum+1) = scalesFirstImPatchIndOffset(tempAlgorithmStepNum);
        if exist('scalePairsFirstIndeces','var')
            scalePairsFirstIndeces(tempAlgorithmStepNum) = scalePairsFirstIndeces(tempAlgorithmStepNum-1);
        else
            scalePairsFirstIndeces(tempAlgorithmStepNum) = 0;
        end
        if tempAlgorithmStepNum<=length(parameters.downscaledSize.multiScale)
            downscaleSize = parameters.downscaledSize.multiScale(tempAlgorithmStepNum);
        else
            noMoreAlgorithmSteps = true;
            break;
        end
    end
    if downscaleSize<returnRelevantImageSize(orgImSize,parameters.downscaledSize.approach)
        algorithmStepNum = tempAlgorithmStepNum;
    end
else
    switch algorithmStepNum
        otherwise
            noMoreAlgorithmSteps = true;
    end
end