%% Parameters that vary accross algorithm steps:
numOfNNs = parameters.NNsearchParams.numOfNNs;
downscaleSize = parameters.downscaledSize.initial;
concatenatedColors = nan(9,3);
configurationsOrderOfMerit = [];    t2DividedByt1_IminusA = [];
downscalingNoiseVar = 0;
airlightRecoveryMethod = 6;
%% Running the multiple-levels algorithm
algorithmStepNum = 1;   scalesFirstImPatchIndOffset = 0;
while 1
    if algorithmStepNum==1 || downscaleSize~=parameters.downscaledSize.initial %only perform this part...
        orgImSize = partialSize(orgImage,1:2);
        %             Avoiding saturated patches. Finding them in the original scale image:
        InvalidDueToPixelsSaturation = double(reshape(max(orgImage==1,[],3),orgImSize));
        if isnan(downscaleSize) ||...
                (parameters.downscaledSize.multiScale(1)>1 && downscaleSize>=returnRelevantImageSize(orgImSize,parameters.downscaledSize.approach)) ||...
                (parameters.downscaledSize.multiScale(1)<=1 && algorithmStepNum==1 &&...
                returnRelevantImageSize(orgImSize,parameters.downscaledSize.approach)<=downscaleSize/sqrt(parameters.downscaledSize.multiScale(1)))
            downscaleSize = returnRelevantImageSize(orgImSize,parameters.downscaledSize.approach);
        else
            InvalidDueToPixelsSaturation = double(imresize2fitDesiredSize(uint8(255*InvalidDueToPixelsSaturation),downscaleSize,parameters.downscaledSize.approach)>0);
            downScaledOrgIm = imresize2fitDesiredSize(uint8(255*orgImage),downscaleSize,parameters.downscaledSize.approach);
            downscalingNoiseVar = var(double(reshape(imresize(downScaledOrgIm,partialSize(orgImage,1:2))-uint8(255*orgImage),[],1))/255);
            orgImage = double(downScaledOrgIm)/255;
        end
        InvalidDueToPixelsSaturation = any(im2col(InvalidDueToPixelsSaturation,parameters.patchSize*ones(1,2),'sliding'));
        if algorithmStepNum==1
            largestScaleUsedImage = orgImage;
        end
        %% Extracting patches
        rescaledOrgIm = orgImage;
        rescaledImDimensions = partialSize(rescaledOrgIm,1:2);
        invalidShoulderPatchesIndices = [];
        tic
        orgImPatches{algorithmStepNum} = rgbim2col(rescaledOrgIm,parameters.patchSize*ones(1,2),'sliding');
        tocAndPrintWhenLarge(2,'Extracting image patches');
        %% Extracting regular patches' STD:
        DCremovedOrgImPatches = bsxfun(@minus,orgImPatches{algorithmStepNum},mean(orgImPatches{algorithmStepNum}));
        executionStruct.functionName = 'permute';
        executionStruct.largeArgumentIndicator = [1,0];
        executionStruct.splittingInputDimension = 2;
        executionStruct.splittedOutputsIndicator = 1;
        executionStruct.splittingOutputDimension = 3;
        [~,DCremovedOrgImPatches] =...
            outOfMemoryResistantExcecution(executionStruct,DCremovedOrgImPatches,[1,3,2]);
        curDir = [];
        executionStruct.functionName = 'std';
        executionStruct.largeArgumentIndicator = 1;
        executionStruct.splittingInputDimension = 2;
        executionStruct.splittedOutputsIndicator = 1;
        executionStruct.splittingOutputDimension = 2;
        [~,orgImPatchesSTDs] = outOfMemoryResistantExcecution(executionStruct,orgImPatches{algorithmStepNum});
        if ~isempty(curDir);  cd(curDir); end;
        %% Parameters calculation:
        patchesMapDimensions = rescaledImDimensions-parameters.patchSize+1; rescaledPatchesMapDimensions(algorithmStepNum,:) = patchesMapDimensions;
        numOfCandidates = min(parameters.NNsearchParams.maxCandidates,size(orgImPatches{algorithmStepNum},2));
        patchPixelsNum = parameters.patchSize^2;
    end
    %% Computing patches' STD:
    tic
    usedPatchesSTD{algorithmStepNum} = orgImPatchesSTDs;
    usedPatchesSTD4Pruning = usedPatchesSTD{algorithmStepNum};
    tocAndPrintWhenLarge(2,'Computing patches'' STD');
    sortedSTDs = sort(usedPatchesSTD{algorithmStepNum},2);
    if isstruct(parameters.NNsearchParams.estimatedNoiseLevelisSTDpercentile)
        [meanSTDhist,histBins] = multiD21Dhist(mean(sortedSTDs,3),1e3);
        [~,maxBin] = max(meanSTDhist);
        meanSTDhist = meanSTDhist(2:end)-meanSTDhist(1:end-1);
        meanSTDhist = conv(meanSTDhist,1/parameters.NNsearchParams.estimatedNoiseLevelisSTDpercentile.diffsAveragingWinSize*...
            ones(1,parameters.NNsearchParams.estimatedNoiseLevelisSTDpercentile.diffsAveragingWinSize),'same');
        histBins = histBins(2:end);
        estimatedNoiseLevelisSTDpercentile(algorithmStepNum) = find(meanSTDhist(maxBin:end)>0,1,'first')+maxBin-1; %#ok<*SAGROW>
        estimatedNoiseLevelisSTDpercentile(algorithmStepNum) = max(parameters.NNsearchParams.estimatedNoiseLevelisSTDpercentile.absoluteMinNoiseSTD,...
            histBins(estimatedNoiseLevelisSTDpercentile(algorithmStepNum)));
        estimatedNoiseLevelisSTDpercentile(algorithmStepNum) = mean(mean(sortedSTDs,3)<=estimatedNoiseLevelisSTDpercentile(algorithmStepNum));
    else
        estimatedNoiseLevelisSTDpercentile(algorithmStepNum) = parameters.NNsearchParams.estimatedNoiseLevelisSTDpercentile;
    end
    estimatedNoiseSTDInLevel{algorithmStepNum} =...
        sortedSTDs(1,round(estimatedNoiseLevelisSTDpercentile(algorithmStepNum)*size(usedPatchesSTD{algorithmStepNum},2)),:);
    estimatedNoiseSTDInLevel{algorithmStepNum} = squeeze(estimatedNoiseSTDInLevel{algorithmStepNum}(1:1:end));
    estimatedNoiseSTDInLevel{algorithmStepNum} = reshape(estimatedNoiseSTDInLevel{algorithmStepNum},1,1,[]);
    %% Pruning by variance
    calculatedValueForPruning = meanOrMaxOrZeroOfValues(usedPatchesSTD4Pruning,3,parameters.NNsearchParams.channelSTD2consider);
    [sortedCalculatedOrgSTD,validCandidatesIndeces] = sort(calculatedValueForPruning,'descend');
    InvalidDueToPixelsSaturation = InvalidDueToPixelsSaturation(validCandidatesIndeces);
    validCandidatesIndeces(InvalidDueToPixelsSaturation) = [];   sortedCalculatedOrgSTD(InvalidDueToPixelsSaturation) = [];
    if isstruct(parameters.pruningThreshold.queriesSTD)
        candidateValidationMap = zeros(1,prod(rescaledPatchesMapDimensions(end,:)));
        candidateValidationMap(validCandidatesIndeces) = 1;
        winSize = round(parameters.pruningThreshold.queriesSTD.winRelativePortion*rescaledPatchesMapDimensions(end,:));
        STDsOfWindows = im2col(reshape(calculatedValueForPruning.*candidateValidationMap,rescaledPatchesMapDimensions(end,:)),winSize,'distinct');
        numOfQueries2Choose = ceil(sum(STDsOfWindows>0)*parameters.pruningThreshold.queriesSTD.upperPercentile);
        [sortedSTDwins,STDwinsOrder] = sort(STDsOfWindows,'descend');
        fixedMinQuerySTD = parameters.pruningThreshold.queriesSTD.noiseLevelWeightRelative2allQueriesPercentile*...
            meanOrMaxOrZeroOfValues(estimatedNoiseSTDInLevel{algorithmStepNum},3,parameters.NNsearchParams.channelSTD2consider)+...
            (1-parameters.pruningThreshold.queriesSTD.noiseLevelWeightRelative2allQueriesPercentile)*...
            percentile(calculatedValueForPruning.',1-parameters.pruningThreshold.queriesSTD.upperPercentile);
        chosenQueryIndicator = false(size(sortedSTDwins));
        chosenQueryIndicator(...
            sortedSTDwins>fixedMinQuerySTD & bsxfun(@le,repmat((1:size(sortedSTDwins,1)).',1,size(sortedSTDwins,2)),numOfQueries2Choose)) = true;
        windowedIndecisMap = im2col(reshape(1:prod(rescaledPatchesMapDimensions(end,:)),rescaledPatchesMapDimensions(end,:)),winSize,'distinct');
        STDwinsOrder = bsxfun(@plus,STDwinsOrder,0:size(STDwinsOrder,1):(numel(STDwinsOrder)-1));
        varPrunedQueryIndeces2keep = windowedIndecisMap(STDwinsOrder(chosenQueryIndicator)).';
    else
        if algorithmStepNum==1
            if parameters.pruningThreshold.queriesSTD>1
                queriesSTDthreshold = percentile(usedPatchesSTD4Pruning(:),parameters.pruningThreshold.queriesSTD/100);
                display(sprintf('Queries STD threshold was set to %.3f',queriesSTDthreshold));
            else
                queriesSTDthreshold = parameters.pruningThreshold.queriesSTD;
            end
        end
        firstHowManyAreQueries = find(sortedCalculatedOrgSTD<=queriesSTDthreshold,1,'first')-1;
        if isempty(firstHowManyAreQueries); firstHowManyAreQueries = length(sortedCalculatedOrgSTD); end;
        firstHowManyAreQueries = max(parameters.NNsearchParams.queriesNumRangeAroundSTD(1),firstHowManyAreQueries);
        varPrunedQueryIndeces2keep = validCandidatesIndeces(1:firstHowManyAreQueries);
    end
    if length(varPrunedQueryIndeces2keep)>parameters.NNsearchParams.queriesNumRangeAroundSTD(2)
        indices2keep = prune2maximizeQueriesPhysicalDistribution(varPrunedQueryIndeces2keep,...
            calculatedValueForPruning(varPrunedQueryIndeces2keep),parameters.NNsearchParams.queriesNumRangeAroundSTD(2),patchesMapDimensions,parameters.patchSize*2+1);
        varPrunedQueryIndeces2keep = varPrunedQueryIndeces2keep(indices2keep);
    end
    if length(varPrunedQueryIndeces2keep)<1
        display(sprintf('Queries number (%d) is below minimum (%d) in algorithm step #%d. Aborting this scale.',...
            length(varPrunedQueryIndeces2keep),parameters.NNsearchParams.minNumOfQueries,algorithmStepNum))
        updateStepParameters
        if noMoreAlgorithmSteps
            validCandidatesIndecesAllScales{algorithmStepNum} = [];
            if all(isnan(concatenatedColors(:)))
                numOfNCpairsCrossingThreshold = []; concatenatedAirlightDists = []; concatenatedColors = [];
            end
            break
        else
            scalesFirstImPatchIndOffset(algorithmStepNum+1) = scalesFirstImPatchIndOffset(algorithmStepNum);
            if algorithmStepNum>1
                if exist('scalePairsFirstIndeces','var')
                    scalePairsFirstIndeces(algorithmStepNum) = scalePairsFirstIndeces(algorithmStepNum-1);
                else
                    scalePairsFirstIndeces(algorithmStepNum) = 0;
                end
            end
            algorithmStepNum = algorithmStepNum+1;
        end
        continue
    end
    validCandidatesIndeces = validCandidatesIndeces(1:min(length(validCandidatesIndeces),numOfCandidates));
    validCandidatesIndeces = sort(validCandidatesIndeces);
    arbitraryPatchIndex4Later = find(validCandidatesIndeces(1:length(varPrunedQueryIndeces2keep))~=varPrunedQueryIndeces2keep,1);
    clear grayOrgImPatches
    NCmap{algorithmStepNum} = calcImageNormalizedCorrelationMap(orgImPatches{algorithmStepNum},rescaledImDimensions,parameters.patchSize);
    preVarNormalizationImPatches{algorithmStepNum} = DCremovedOrgImPatches;
    preVarNormalizationImPatches4search = preVarNormalizationImPatches{algorithmStepNum};
    %     if parameters.NNsearchParams.colorCoupledNormalization
    NormalizedDBPatches{algorithmStepNum} = reshape(preVarNormalizationImPatches4search(:,:,validCandidatesIndeces),...
        size(preVarNormalizationImPatches4search,2)*parameters.patchSize^2,[]);
    NormalizedQueryPatches{algorithmStepNum} = reshape(preVarNormalizationImPatches4search(:,:,varPrunedQueryIndeces2keep),...
        size(preVarNormalizationImPatches4search,2)*parameters.patchSize^2,[]);
    
    executionStruct.functionName = 'std';
    executionStruct.largeArgumentIndicator = 1;
    executionStruct.splittingInputDimension = 2;
    executionStruct.splittedOutputsIndicator = 1;
    executionStruct.splittingOutputDimension = 2;
    [~,curAlgStepNormalizedDBPatchesSTD] =...
        outOfMemoryResistantExcecution(executionStruct,NormalizedDBPatches{algorithmStepNum});
    NormalizedDBPatches{algorithmStepNum} = bsxfun(@rdivide,NormalizedDBPatches{algorithmStepNum},...
        bsxfun(@max,curAlgStepNormalizedDBPatchesSTD,mean(estimatedNoiseSTDInLevel{algorithmStepNum})));
    NormalizedQueryPatches{algorithmStepNum} = bsxfun(@rdivide,NormalizedQueryPatches{algorithmStepNum},...
        bsxfun(@max,std(NormalizedQueryPatches{algorithmStepNum}),mean(estimatedNoiseSTDInLevel{algorithmStepNum})));
    croppedOrgIm = patchLikeCropAnImage(rescaledOrgIm,parameters.patchSize);
    if algorithmStepNum==1
        DBFoldIndices = [];
        [curScaleDBFoldIndices,DBpatchesFoldAssignment{1},initialTransmissionGuess,foldsAssignmentMap] = assignPatches2Folds(validCandidatesIndeces,...
            patchesMapDimensions,parameters.NNsearchParams.numOfNNsearchFolds,'brightnessVarCombination',squeeze(mean(usedPatchesSTD{algorithmStepNum},3)),...
            rescaledOrgIm,[],parameters.patchSize);
    else
        [curScaleDBFoldIndices,DBpatchesFoldAssignment{algorithmStepNum}] = assignPatches2Folds(validCandidatesIndeces,patchesMapDimensions,...
            parameters.NNsearchParams.numOfNNsearchFolds,'brightnessVarCombination',squeeze(mean(usedPatchesSTD{algorithmStepNum},3)),rescaledOrgIm,[],parameters.patchSize,round(imresize(foldsAssignmentMap,patchesMapDimensions)));
    end
    DBFoldIndices = [DBFoldIndices,curScaleDBFoldIndices]; %#ok<*AGROW>
    [~,~,queryFoldIndeces{algorithmStepNum}] = intersect(varPrunedQueryIndeces2keep,validCandidatesIndeces,'stable');
    queryFoldIndeces{algorithmStepNum} = DBpatchesFoldAssignment{algorithmStepNum}(queryFoldIndeces{algorithmStepNum});
    varPrunedQueryIndeces2keepAllScales{algorithmStepNum} = varPrunedQueryIndeces2keep+scalesFirstImPatchIndOffset(algorithmStepNum);
    validCandidatesIndecesAllScales{algorithmStepNum} = validCandidatesIndeces+scalesFirstImPatchIndOffset(algorithmStepNum);
    scalesFirstImPatchIndOffset(algorithmStepNum+1) = scalesFirstImPatchIndOffset(algorithmStepNum)+size(preVarNormalizationImPatches{algorithmStepNum},3);
    updateStepParameters
    if noMoreAlgorithmSteps
        break
    else
        algorithmStepNum = algorithmStepNum+1;
    end
end
numOfScales = algorithmStepNum;
scalesFirstQueryIndOffset = cumsum([0,cellfun('length',queryFoldIndeces)]); %scalesFirstImPatchIndOffset = [0,cellfun('size',preVarNormalizationImPatches,3)];
nonEmptyScalesIndicator = ~cellfun('isempty',NormalizedDBPatches);
reshapedpreVarNormalizationImPatches =...
    reshape(cell2mat(reshape(preVarNormalizationImPatches,1,1,[])),prod(partialSize(preVarNormalizationImPatches{1},1:2)),scalesFirstImPatchIndOffset(end));
scalesFirstImPatchIndOffset = [0,cumsum(cellfun('size',preVarNormalizationImPatches,3).*~cellfun('isempty',preVarNormalizationImPatches))];
clear preVarNormalizationImPatches
%% Finding NNs:
display('Searching patches for nearest neighbors...')
KDsearchTotalTime = 0;
numNN2search4 = parameters.NNsearchParams.maxNNsFromEachFold;
maxNumOfFoldsPerScale = parameters.NNsearchParams.numOfNNsearchFolds-2;
chosenIndeces = nan(numOfScales*numNN2search4*maxNumOfFoldsPerScale,sum(cellfun('length',queryFoldIndeces)));
if (parameters.NNsearchParams.numOfNNsearchFolds-3)*parameters.NNsearchParams.maxNNsFromEachFold<numOfNNs
    error('Too many desired NNs Or too little NNs per fold Or too little folds Or don''t restrict NNs to ''far''')
end
clear bestMatches complianceWithAirlight
allScalesDB2nonEmptyScalesDBmapping = nan(size(nonEmptyScalesIndicator));   allScalesDB2nonEmptyScalesDBmapping(nonEmptyScalesIndicator) = 1:sum(nonEmptyScalesIndicator);
rng(1);
for queryScaleNum = 1:length(NormalizedQueryPatches)
    if isempty(NormalizedQueryPatches{queryScaleNum});  continue;   end;
    for queryFoldNum = 1:parameters.NNsearchParams.numOfNNsearchFolds
        if ~any(queryFoldIndeces{queryScaleNum}==queryFoldNum)
            continue;
        end
        curFoldQueryInd = find(queryFoldIndeces{queryScaleNum}==queryFoldNum);
        DBfoldNumsPerQueryFold = [1:(queryFoldNum-2) (queryFoldNum+2):size(DBFoldIndices,1)];
        perQueryDuration = 0;
        maxDBscale2usedAsDB = queryScaleNum;
        for DBscaleNum = 1:maxDBscale2usedAsDB
            if isempty(NormalizedDBPatches{DBscaleNum});    continue;   end;
            DBfoldCounterInsideScale = 0;
            for DBfoldNum = DBfoldNumsPerQueryFold
                DBfoldCounterInsideScale = DBfoldCounterInsideScale+1;
                curFoldDBInd = DBFoldIndices{DBfoldNum,allScalesDB2nonEmptyScalesDBmapping(DBscaleNum)};
                tempRowIndeces = (DBscaleNum-1)*maxNumOfFoldsPerScale*parameters.NNsearchParams.maxNNsFromEachFold+...
                    ((DBfoldCounterInsideScale-1)*parameters.NNsearchParams.maxNNsFromEachFold+1:DBfoldCounterInsideScale*parameters.NNsearchParams.maxNNsFromEachFold);
                %                             end
                if length(curFoldDBInd)<numNN2search4
                    continue;
                end
                searchEngine = 'ANN';%'ANN','KDtree'
                currentGlobalQueryIndeces = curFoldQueryInd+scalesFirstQueryIndOffset(queryScaleNum);
                switch searchEngine
                    case 'ANN'
                        tic
                        if ~exist('kdtree','var') || any(size(kdtree)<[DBscaleNum,DBfoldNum]) || isempty(kdtree{DBscaleNum,DBfoldNum})
                            kdtree{DBscaleNum,DBfoldNum} = ann(single(NormalizedDBPatches{DBscaleNum}(:,curFoldDBInd)));
                        end
                        chosenIndeces(tempRowIndeces,currentGlobalQueryIndeces) =...
                            ksearch(kdtree{DBscaleNum,DBfoldNum},single(NormalizedQueryPatches{queryScaleNum}(:,curFoldQueryInd)),numNN2search4,4,0);
                        searchDuration = toc;
                    case 'KDtree'
                        executionStruct.functionName = 'complex_KDSearch';
                        executionStruct.largeArgumentIndicator = [false,true,false,false];
                        executionStruct.splittingInputDimension = 2;
                        executionStruct.splittedOutputsIndicator = [true,false,false,true];
                        executionStruct.splittingOutputDimension = 2;
                        [~,chosenIndeces(tempRowIndeces,currentGlobalQueryIndeces),~,searchDuration] =...
                            outOfMemoryResistantExcecution(executionStruct,single(NormalizedDBPatches{DBscaleNum}(:,curFoldDBInd)),...
                            single(NormalizedQueryPatches{queryScaleNum}(:,curFoldQueryInd)),numNN2search4,1);
                end
                perQueryDuration = perQueryDuration+searchDuration;
                chosenIndeces(tempRowIndeces,currentGlobalQueryIndeces) = curFoldDBInd(chosenIndeces(tempRowIndeces,currentGlobalQueryIndeces));
            end
        end
        KDsearchTotalTime = KDsearchTotalTime+perQueryDuration;
    end
end
if strcmp(searchEngine,'ANN')
    for treeNum1 = 1:size(kdtree,1)
        for treeNum2 = 1:size(kdtree,2)
            temp = close(kdtree{treeNum1,treeNum2});
        end
    end
end
clear NormalizedDBPatches
tocAndPrintWhenLarge(2,'Searched for nearest neighbors',KDsearchTotalTime);
%% Picking best NNs:
chosenIndecesOffsetVect = cumsum(cellfun('size',validCandidatesIndecesAllScales,2));
validCandidatesIndecesAllScales = cell2mat(validCandidatesIndecesAllScales);    varPrunedQueryIndeces2keepAllScales = cell2mat(varPrunedQueryIndeces2keepAllScales);
chosenIndecesOffsetVect = reshape(repmat([0,chosenIndecesOffsetVect(1:end-1)],maxNumOfFoldsPerScale*parameters.NNsearchParams.maxNNsFromEachFold,1),[],1);
chosenIndeces = bsxfun(@plus,chosenIndeces,chosenIndecesOffsetVect);
chosenIndecesNansLocations = isnan(chosenIndeces);
chosenIndeces(chosenIndecesNansLocations) = arbitraryPatchIndex4Later;%Getting rid of nans (in case there is an unequal number of NNs for each query, so the remaining are nans) by putting there the arbitrary NN 1.
queriesCorrespondingMinSTDs = zeros(scalesFirstQueryIndOffset(end),3);
queriesCorrespondingMinSTDs(1+scalesFirstQueryIndOffset([true,nonEmptyScalesIndicator(1:end-1)]),:) =...
    cat(2,estimatedNoiseSTDInLevel{1},diff(cell2mat(estimatedNoiseSTDInLevel),1,2));
queriesCorrespondingMinSTDs = shiftdim(cumsum(queriesCorrespondingMinSTDs),-1);
DBcandsCorrespondingMinSTD = zeros(1,numOfScales,3);    DBcandsCorrespondingMinSTD(1,nonEmptyScalesIndicator,:) = cell2mat(estimatedNoiseSTDInLevel);
separatelyLimitEachSTDinNC = true;
minDBSTDs = reshape(repmat(DBcandsCorrespondingMinSTD,maxNumOfFoldsPerScale*...
    parameters.NNsearchParams.maxNNsFromEachFold,1),[],1,3);
minQuerySTDs = queriesCorrespondingMinSTDs;
executionStruct.functionName = 'calcNormalizedCorrelation';
executionStruct.largeArgumentIndicator = [false,true,true,false,false,true,false];
executionStruct.splittingInputDimension = 2;
executionStruct.splittedOutputsIndicator = [true,true];
executionStruct.splittingOutputDimension = 2;
executionStruct.splittingFactor = ceil(numel(chosenIndeces)/1e6)*10.^(0:3);
[~,curStageNNnormalizedCorrelationess,targetPatchesVar] = outOfMemoryResistantExcecution(executionStruct,...
    reshapedpreVarNormalizationImPatches,reshape(validCandidatesIndecesAllScales(chosenIndeces),size(chosenIndeces)),...
    varPrunedQueryIndeces2keepAllScales,parameters.patchSize^2,minDBSTDs,minQuerySTDs,separatelyLimitEachSTDinNC);
curStageNNnormalizedCorrelationess(chosenIndecesNansLocations) = -1;% Making sure the random patches added 4 lines ago don't get chosen, because I don't have curStageNNdistances for them. Putting 0 here...
clear minDBSTDs minQuerySTDs
tocAndPrintWhenLarge(2,'Calculating neighbors'' normalized correlations');
multiFoldEncouragingFactors = parameters.NNsearchParams.multiFoldEncouragingFactor.^(0:numNN2search4-1);
multiFoldEncouragingFactors = repmat(multiFoldEncouragingFactors',maxNumOfFoldsPerScale*numOfScales,1);
[~,NNsWeightedNCOrder] = sort(bsxfun(@times,curStageNNnormalizedCorrelationess,multiFoldEncouragingFactors),'descend');% Sorting according to norm correlation weighted by multi-folds encouragement
if ~isempty(parameters.pruningThreshold.minNormalizedCorrelation)
    targetPatchesVar(curStageNNnormalizedCorrelationess<parameters.pruningThreshold.minNormalizedCorrelation) =...
        targetPatchesVar(curStageNNnormalizedCorrelationess<parameters.pruningThreshold.minNormalizedCorrelation)+max(targetPatchesVar(:));
end
[~,targetPatchesVarOrder] = sort(targetPatchesVar);
executionStruct.functionName = 'mergeMultipleSorts';
executionStruct.largeArgumentIndicator = true;
executionStruct.splittingInputDimension = 3;
executionStruct.splittedOutputsIndicator = true;
executionStruct.splittingOutputDimension = 1;
[~,NNpreferenceOrder] = outOfMemoryResistantExcecution(executionStruct,...
    cat(1,shiftdim(NNsWeightedNCOrder,-1),shiftdim(targetPatchesVarOrder,-1)));
clear NNsWeightedNCOrder targetPatchesVar targetPatchesVarOrder
NNpreferenceOrder = NNpreferenceOrder.';
NNpreferenceOrder = NNpreferenceOrder(1:numOfNNs-1,:);% Keeping only the desired number of neighbors
tempIndeces = sub2ind(size(chosenIndeces),NNpreferenceOrder,repmat(1:size(chosenIndeces,2),numOfNNs-1,1));% computing global indeces to keep
chosenIndeces = chosenIndeces(tempIndeces);%curStageNNdistances = curStageNNdistances(tempIndeces);
curStageNNnormalizedCorrelationess = curStageNNnormalizedCorrelationess(tempIndeces);% keeping only the chosen indeces
[curStageNNnormalizedCorrelationess,NNpreferenceOrder] = sort(curStageNNnormalizedCorrelationess,'descend');% Rearranging the columns to be in descnding (non-weighted) normalized correlation order
tempIndeces = sub2ind(size(chosenIndeces),NNpreferenceOrder,repmat(1:size(chosenIndeces,2),numOfNNs-1,1));
chosenIndeces = chosenIndeces(tempIndeces);%curStageNNdistances = curStageNNdistances(tempIndeces);
clear bestMatches
NNnormalizedCorrelation = [ones(1,size(curStageNNnormalizedCorrelationess,2));curStageNNnormalizedCorrelationess];
bestMatches(2:numOfNNs,:) = reshape(validCandidatesIndecesAllScales(chosenIndeces),size(chosenIndeces));
bestMatches(1,:) = varPrunedQueryIndeces2keepAllScales;
[~,scaleNum] = max(bsxfun(@gt,bestMatches(:),scalesFirstImPatchIndOffset(1:end-1)) & bsxfun(@le,bestMatches(:),scalesFirstImPatchIndOffset(2:end)),[],2);
bestMatchesCorrespondingScales = reshape(scaleNum,size(bestMatches));
clear NormalizedDBPatches
patchGroupNum = repmat(1:size(bestMatches,2),size(bestMatches,1),1);
bestMatches = cat(1,repmat(bestMatches(1,:),1,size(bestMatches,1)-1),reshape(bestMatches(2:end,:).',1,[]));
bestMatchesCorrespondingScales = cat(1,repmat(bestMatchesCorrespondingScales(1,:),1,size(bestMatchesCorrespondingScales,1)-1),...
    reshape(bestMatchesCorrespondingScales(2:end,:).',1,[]));
patchGroupNum = cat(1,repmat(patchGroupNum(1,:),1,size(patchGroupNum,1)-1),reshape(patchGroupNum(2:end,:).',1,[]));
NNnormalizedCorrelation = reshape(NNnormalizedCorrelation(2:end,:).',1,[]);

executionStruct.functionName = 'reshape4outOfMemoryResistantExcecution';
executionStruct.largeArgumentIndicator = [1,0,0];
executionStruct.splittingInputDimension = 2;
executionStruct.splittedOutputsIndicator = 1;
executionStruct.splittingOutputDimension = 4;
maxSplittingParts = size(bestMatches,2)/gcd(size(bestMatches,2),2^ceil(log2(size(bestMatches,2))));
executionStruct.splittingFactor = divisor(maxSplittingParts);
executionStruct.splittingFactor = executionStruct.splittingFactor(1:min(4,length(executionStruct.splittingFactor)));
[~,DClessOrAlessChosenPatches] =...
    outOfMemoryResistantExcecution(executionStruct,reshapedpreVarNormalizationImPatches(:,bestMatches),...
    [1,1,2,2],[parameters.patchSize^2,3,2,nan]);

clear reshapedpreVarNormalizationImPatches
executionStruct.functionName = 'permute';
executionStruct.largeArgumentIndicator = [1,0];
executionStruct.splittingInputDimension = 4;
executionStruct.splittedOutputsIndicator = 1;
executionStruct.splittingOutputDimension = 3;
[~,DClessOrAlessChosenPatches] =...
    outOfMemoryResistantExcecution(executionStruct,DClessOrAlessChosenPatches,[1,3,4,2]);

initialTratioCalculationMethod = ternaryOperator(strcmp(parameters.tRatioCalculation,'factorization'),'norm2Ratios',parameters.tRatioCalculation);
t2DividedByt1 = reshape(computePatchesRatio(DClessOrAlessChosenPatches,initialTratioCalculationMethod),size(bestMatches)-[1,0]);
bestMatchesIndeces2keep = find(t2DividedByt1<1); %this will also prevent the case of dual pair appearance, except when t-ratios are not calculated in pairs and when pairs with...
%     NC below threshold manged to enter the NNs list (because there were
%     no better ones).
%         Discarding pairs that appear twice (in opposite order):
[~,bestMatchesIndeces2keep_temp] = unique(sort(bestMatches(:,bestMatchesIndeces2keep)).','rows','stable');
bestMatchesIndeces2keep = bestMatchesIndeces2keep(bestMatchesIndeces2keep_temp);
DClessOrAlessChosenPatches = DClessOrAlessChosenPatches(:,:,bestMatchesIndeces2keep,:); %t2DividedByt1 = t2DividedByt1(bestMatchesIndeces2keep);
bestMatches = bestMatches(:,bestMatchesIndeces2keep);   NNnormalizedCorrelation = NNnormalizedCorrelation(:,bestMatchesIndeces2keep);
patchGroupNum = patchGroupNum(1,bestMatchesIndeces2keep);   bestMatchesCorrespondingScales = bestMatchesCorrespondingScales(:,bestMatchesIndeces2keep);
t2DividedByt1 = reshape(computePatchesRatio(...
    DClessOrAlessChosenPatches,parameters.tRatioCalculation,patchGroupNum,NNnormalizedCorrelation,parameters.pruningThreshold.minNormalizedCorrelation),size(bestMatches)-[1,0]);
bestMatches_perScalePatchIndeces = bestMatches-scalesFirstImPatchIndOffset(bestMatchesCorrespondingScales);
if ~isempty(NCmap)
    bestMatchesCorrespondingShiftsNC = nan(size(bestMatchesCorrespondingScales));
    for scaleNum = 1:numOfScales
        if any(bestMatchesCorrespondingScales(:)==scaleNum)
            bestMatchesCorrespondingShiftsNC(bestMatchesCorrespondingScales==scaleNum) = NCmap{scaleNum}(bestMatches_perScalePatchIndeces(bestMatchesCorrespondingScales==scaleNum));
        end
    end
else
    bestMatchesCorrespondingShiftsNC = [];
end
clear NormalizedQueryPatches
clear DCremovedOrgImPatches preVarNormalizationImPatches4search
orgImPatches = cell2mat(orgImPatches);
originalStateChosenPatches = reshape(orgImPatches(:,bestMatches,:),[size(orgImPatches,1),size(bestMatches),size(orgImPatches,3)]);
clear orgImPatches
%% Thresholds & logic based filtering
thresholdBasedValidity = all(t2DividedByt1>0,3);
if ~isempty(parameters.pruningThreshold.minNormalizedCorrelation)% Minimum normalized correlation
    thresholdBasedValidity = thresholdBasedValidity & NNnormalizedCorrelation>=returnLowerThreshold(NNnormalizedCorrelation(thresholdBasedValidity),...
        parameters.pairsPruning.globalRate(1),parameters.pruningThreshold.minNormalizedCorrelation);
end
if ~isempty(parameters.pairsPruning.allowed_t2Tot1_ratioRange)% t-ratios range:
    thresholdBasedValidity = thresholdBasedValidity &...
        all(abs(log(t2DividedByt1))>parameters.pairsPruning.allowed_t2Tot1_ratioRange(1) & abs(log(t2DividedByt1))<parameters.pairsPruning.allowed_t2Tot1_ratioRange(2),3);
end
betasRatioWiseValidity = nan(size(NNnormalizedCorrelation));
NormalizedCorrWiseValidity = NNnormalizedCorrelation;
aprioriInvalidity = false(size(betasRatioWiseValidity));

max_airlightEstimationIterations = 3;
minPairsNumberShrinkFactor = .5;

overallNNsValidity = true(size(NormalizedCorrWiseValidity));    AirlightEstimationIterNum = 0;  alreadyPrunedByRecScore = false;
while AirlightEstimationIterNum<max_airlightEstimationIterations
    AirlightEstimationIterNum = AirlightEstimationIterNum+1;    previousOverallNNsValidity = overallNNsValidity;
    [~,overallNNsValidityVals] = computeOverallValidity(parameters.pairsPruning,parameters.pairsPruning.globalRate(1),...
        betasRatioWiseValidity,NormalizedCorrWiseValidity,reshape(t2DividedByt1,partialSize(t2DividedByt1,2:3)),aprioriInvalidity);
    overallNNsValidityVals(~thresholdBasedValidity) = nan;
    if AirlightEstimationIterNum>1
        airlightRecoveryMethod = tempHoldAirlightRecoveryMethod;
        varOfPCAofIndAsComponents = eig(cov(individualPairAirlight));
        numOfPairsUsed = parameters.pairsPruning.globalRate;
        if exist('tempHoldPrioritizingMethod','var')
            parameters.pairsPruning.overlapPruningPrioritizingMethod = tempHoldPrioritizingMethod;
        end
        pruningVals.reoccurrenceVals = nan(1,length(overallNNsValidityVals));
        pruningVals.reoccurrenceVals(overallNNsValidity) = reoccurrenceScore{parameters.airlightRecovery.reoccurrenceScore.methodNum};
        NNsValidityValsOfLatestValidPairs = overallNNsValidityVals(latestIterValidPairs);
        overallNNsValidityVals = nan(size(overallNNsValidityVals)); overallNNsValidityVals(latestIterValidPairs) =...
            NNsValidityValsOfLatestValidPairs;
        upperReccurrenceScoreThreshold = -1*returnLowerThreshold(-reoccurrenceScore{parameters.airlightRecovery.reoccurrenceScore.methodNum},...
            parameters.pairsPruning.globalRate(1),-parameters.airlightRecovery.reoccurrenceScore.allowedRange(2));
        latestValidPairs2discard = reoccurrenceScore{parameters.airlightRecovery.reoccurrenceScore.methodNum}<...
            parameters.airlightRecovery.reoccurrenceScore.allowedRange(1) |...
            reoccurrenceScore{parameters.airlightRecovery.reoccurrenceScore.methodNum}>upperReccurrenceScoreThreshold;
        if all(latestValidPairs2discard);   error('Recurrence score discarded all pairs');  end;
        overallNNsValidityVals(latestIterValidPairs(latestValidPairs2discard)) = nan;
        individualPairAirlight(reoccurrenceScore{parameters.airlightRecovery.reoccurrenceScore.methodNum}<parameters.airlightRecovery.reoccurrenceScore.allowedRange(1) |...
            reoccurrenceScore{parameters.airlightRecovery.reoccurrenceScore.methodNum}>parameters.airlightRecovery.reoccurrenceScore.allowedRange(2),:) = [];
%         overallNNsValidity(isnan(overallNNsValidityVals)) = false;
%         latestIterValidPairs = find(overallNNsValidity);
        alreadyPrunedByRecScore = true;
    else
        if strcmp(parameters.pairsPruning.overlapPruningPrioritizingMethod,'reoccurrenceScore') ||...
                strcmp(parameters.pairsPruning.overlapPruningPrioritizingMethod,'reoccurrence_complexity_meansVariety_tRatios')
            tempHoldPrioritizingMethod = parameters.pairsPruning.overlapPruningPrioritizingMethod;
            parameters.pairsPruning.overlapPruningPrioritizingMethod = 'complexity_meansVariety_tRatios';
        end
        tempHoldAirlightRecoveryMethod = airlightRecoveryMethod;
        if airlightRecoveryMethod==6 %|| airlightRecoveryMethod==9%The noise model calculation of A cannot handle a large number of pairs, so I'm using the no-noise model for the initial outlier rejection.
            airlightRecoveryMethod = 4;
        end
        t2DividedByt1EstimationDiffs = -overallNNsValidityVals;
        numOfPairsUsed = parameters.pairsPruning.globalRate(1);
        numOfPairsUsed(2) = 3e4;
    end
    %% Assiging validity indicator for this A estimation iteration
    pruningVals.overallNNsValidityVals = overallNNsValidityVals;    pruningVals.t2DividedByt1EstimationDiffs = t2DividedByt1EstimationDiffs;
    pruningVals.bestMatchesCorrespondingShiftsNC = bestMatchesCorrespondingShiftsNC;    pruningVals.tRatios = mean(t2DividedByt1,3);
    overallNNsValidity = prioritizeNoOverlapAndPrunePairs(pruningVals,originalStateChosenPatches(:,1,:,:),bestMatches_perScalePatchIndeces,bestMatchesCorrespondingScales,...
        rescaledPatchesMapDimensions,parameters.pairsPruning.overlapPruningPrioritizingMethod,parameters.airlightRecovery.maxAllowedPairsOverlap,parameters.patchSize,...
        numOfPairsUsed,parameters.pruningThreshold.minNormalizedCorrelation,AirlightEstimationIterNum>1); %Not pruning overlapping patches in the first iteration of airlight...
    %         estimation, because it is only used for individual As, and the result of their clustering process is used for the t map propagation. Only after...
    % the clustering I prune the overlapping ones for the global A estimation
    latestIterValidPairs = find(overallNNsValidity);
    if ~any(overallNNsValidity)
        if noMoreAlgorithmSteps
            if all(isnan(concatenatedColors(:)))
                numOfNCpairsCrossingThreshold = []; concatenatedAirlightDists = []; concatenatedColors = [];
            end
        else
            if algorithmStepNum>1
                if exist('scalePairsFirstIndeces','var')
                    scalePairsFirstIndeces(algorithmStepNum) = scalePairsFirstIndeces(algorithmStepNum-1);
                else
                    scalePairsFirstIndeces(algorithmStepNum) = 0;
                end
            end
        end
        break
    end
    %% Airlight estimation
    %% Calculating estimated t1 values:
    if any(airlightRecoveryMethod==6 | airlightRecoveryMethod==8)
        bestMatchesInHighestScale = rescalePatchMapLocations2HighestLevel(bestMatches_perScalePatchIndeces(:,overallNNsValidity),...
            bestMatchesCorrespondingScales(:,overallNNsValidity),rescaledPatchesMapDimensions,parameters.patchSize);
        patchesTestimationOrNC = [];
    else
        patchesTestimationOrNC = [];
    end
    %% Estimating A:
    overallNNsValidity_firstEstimationOfThisStep = overallNNsValidity;
    if isempty(parameters.downscaledSize.multiScale)
        scalePairsFirstIndeces = [1,length(NNnormalizedCorrelation(overallNNsValidity))+1];
    end
    [curAirlightColors,individualPairAirlight,pairsWeights_validPatches] =...
        estimateAirlightColor_OuterFunc(DClessOrAlessChosenPatches,originalStateChosenPatches,overallNNsValidity,...
        patchesTestimationOrNC,t2DividedByt1,airlightRecoveryMethod);
    individualPairAirlight = individualPairAirlight(:,:,airlightRecoveryMethod);
    concatenatedColors(1:size(curAirlightColors,1),:) = curAirlightColors;
    reoccurrenceScore = calcAreoccurrenceScore(originalStateChosenPatches(:,:,overallNNsValidity,:),...
        individualPairAirlight,parameters.airlightRecovery.reoccurrenceScore.omega4lowerBound);
    reoccurrenceScore{4}(any(individualPairAirlight<0 | individualPairAirlight>1,2)) = max(reoccurrenceScore{4})+1;
    %% Evaluating Airlight estimation:
    curAirlightVectEst = concatenatedColors(airlightRecoveryMethod,:);
    airlight2complyWith = curAirlightVectEst;
    max_airlightEstimationIterations = ternaryOperator((~isempty(parameters.airlightRecovery.reoccurrenceScore.methodNum)),2,AirlightEstimationIterNum);
    numOfNCpairsCrossingThreshold = sum(NNnormalizedCorrelation(overallNNsValidity)>parameters.pruningThreshold.minNormalizedCorrelation);
    bestMatches_validPatches = bestMatches_perScalePatchIndeces(:,overallNNsValidity);  % rescaledPatchesMapDimensions(algorithmStepNum,:) = patchesMapDimensions;
    bestMatchesCorrespondingScales_validPatches = bestMatchesCorrespondingScales(:,overallNNsValidity);
    t2DividedByt1_ValidPatches = t2DividedByt1(:,overallNNsValidity,:);
    NNnormalizedCorrelation_validPatches = NNnormalizedCorrelation(overallNNsValidity);
end
concatenatedColors = concatenatedColors(airlightRecoveryMethod,:);
updatedCalculatedAirlight = concatenatedColors;