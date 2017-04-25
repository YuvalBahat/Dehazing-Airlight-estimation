function    [overallValidity_Binary,overallValidity_Real] = computeOverallValidity(parameters,globalSurvivalRate,betasRatio,NormalizedCorr,t2DividedByt1,aprioriInvalidity)
overallValidity_Real = ones(size(NormalizedCorr));
overallValidity_Real(aprioriInvalidity) = 0;
% if ~isempty(parameters.('withinBetaQuotientRangeProbability'))
%     betasRatio = histeq(betasRatio/max(betasRatio(:)),1e5);
%     overallValidity_Real = overallValidity_Real.*(betasRatio.^(1/parameters.('withinBetaQuotientRangeProbability')));
% end
% switch parameters.pruningMethod
%     case 'NC_betaRatios'
%         overallValidity_Real = overallValidity_Real.*(betasRatio)/max(betasRatio(:)).*(NormalizedCorr/max(NormalizedCorr(:)).*(NormalizedCorr>=0)).^10;
%     case 'betaRatios'
% %         betasRatio = histeq(betasRatio/max(betasRatio(:)),1e5);
%         overallValidity_Real = overallValidity_Real.*(betasRatio);
%     case 'NC'
        overallValidity_Real = overallValidity_Real & (NormalizedCorr>=0);
        NormalizedCorr(NormalizedCorr<0) = 0;
%         NormalizedCorr = NormalizedCorr/max(NormalizedCorr(:));
        overallValidity_Real = overallValidity_Real.*(NormalizedCorr);
%     case 'minT2DividedByT1TimesNC'
% %         overallValidity_Real = bsxfun(@times,sigmf(min((1-t2DividedByt1(:,:,:)),[],2),[8,.5]),NormalizedCorr.').';
% %         overallValidity_Real = bsxfun(@times,min((1-t2DividedByt1(:,:,:)),[],2),NormalizedCorr.').';
%         overallValidity_Real = bsxfun(@times,bsxfun(@min,min((1-t2DividedByt1(:,:,:)),[],2),.3),NormalizedCorr.').';
% end
% if ~isempty(parameters.('NNnormCorrelation'))
%     overallValidity_Real = overallValidity_Real & (NormalizedCorr>=0);
%     NormalizedCorr(NormalizedCorr<0) = 0;
%     NormalizedCorr = NormalizedCorr/max(NormalizedCorr(:));
%     overallValidity_Real = overallValidity_Real.*(NormalizedCorr.^(1/parameters.('NNnormCorrelation')));
% end
% % if ~isempty(parameters.('normsQuotientRange'))
% %     normRatios = normRatios/max(normRatios(:));
% %     overallValidity_Real = overallValidity_Real.*(normRatios.^(1/parameters.('normsQuotientRange')));
% % end
% if ~isempty(parameters.('channelsBetaRatios_1NN'))
%     betasRatio = histeq(betasRatio/max(betasRatio(:)),1e5);
%     overallValidity_Real = overallValidity_Real.*(betasRatio.^(1/parameters.('channelsBetaRatios_1NN')));
% end
if globalSurvivalRate<=1 && ~globalSurvivalRate==0% Setting the valid portion:
    validityThreshold = returnThresholdValue(overallValidity_Real(~aprioriInvalidity(:)),globalSurvivalRate,'lower');
    %     notchBelowThreshold = validityThreshold(round(globalSurvivalRate*length(validityThreshold)*mean(~aprioriInvalidity(:)))+1);
    %     validityThreshold = validityThreshold(end:-1:1);
    %     thresholdingIndex = find(validityThreshold>notchBelowThreshold,1,'first');
    %     validityThreshold = validityThreshold(thresholdingIndex);
elseif globalSurvivalRate==0
    validityThreshold = min(overallValidity_Real(:));
else %Setting the number of valid NNs:
    validityThreshold = sort(overallValidity_Real(:),'descend');
    validityThreshold = validityThreshold(min(globalSurvivalRate,length(validityThreshold)));
end
overallValidity_Binary = overallValidity_Real>=validityThreshold;