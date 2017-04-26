function [indexInMatrix,orderOfIdentifiers2Use,maxIdentifierRecurrence,numOfUniqueIdentifiers] =...
    returnGatherByCommonIdentifiersStructer(identifierNum2UniqeIdentifierIndex,numOfUniqueIdentifiers,activeIdentifiersOnly)
% a 2xN array identifierNum2UniqeIdentifierIndex of mapping between
% identifiers to their locations in a unique version, is transformed so that
% identifiers mapped to the same uinque location can be summed.
% Inputs:
% identifierNum2UniqeIdentifierIndex - MxN, mapping between identifiers to
% their index in a unique locations vector. M doesn't really matter. NOTE
% that orderOfIdentifiers2Use works first on the 2nd dimension of
% identifierNum2UniqeIdentifierIndex.
% numOfUniqueIdentifiers - Length of the unique locations vector
% activeIdentifiersOnly - true - In case identifierNum2UniqeIdentifierIndex does
% not include all integers in the range
% 1:max(identifierNum2UniqeIdentifierIndex), the matrix of unique
% summations does not have numOfUniqueIdentifiers columns, but instead has
% length(unique(identifierNum2UniqeIdentifierIndex)) columns. The indeces
% in indexInMatrix are chnaged accordingly.

% Outputs:
% The matrix for unique identifiers summation is of size maxIdentifierRecurrence x
% numOfUniqueIdentifiers. indexInMatrix is the index in the matrix of
% identifierNum2UniqeIdentifierIndex, ordered according to
% orderOfIdentifiers2Use, where the order is not the regular one, but
% rather works on the first and then the second row of
% identifierNum2UniqeIdentifierIndex.
% The function will ignore nans in identifierNum2UniqeIdentifierIndex so
% that orderOfIdentifiers2Use addresses only the non-nan indices of identifierNum2UniqeIdentifierIndex
% If activeIdentifiersOnly, the columns of the resulting matrix keep the original order of unique(identifierNum2UniqeIdentifierIndex,'stable').
% If ~activeIdentifiersOnly, then the columns are in an ascending order of the
% values in identifierNum2UniqeIdentifierIndex.
if nargin<3
    activeIdentifiersOnly = false;
end
nanLocationsInInputVector = isnan(identifierNum2UniqeIdentifierIndex.');
if activeIdentifiersOnly
    [unique2identiferLocationMapping,~,identifierNum2UniqeIdentifierIndex] = unique(reshape(identifierNum2UniqeIdentifierIndex.',[],1),'stable');
    numOfUniqueIdentifiers = length(unique2identiferLocationMapping(~isnan(unique2identiferLocationMapping)));
    identifierNum2UniqeIdentifierIndex(nanLocationsInInputVector) = nan;
else
    identifierNum2UniqeIdentifierIndex = reshape(identifierNum2UniqeIdentifierIndex.',[],1);
end
[orderedIdenifiersLocations,orderOfIdentifiers2Use] = sort(identifierNum2UniqeIdentifierIndex);
orderedIdenifiersLocations(isnan(orderedIdenifiersLocations)) = [];
if activeIdentifiersOnly
    [~,~,orderedIdenifiersLocations] = unique(orderedIdenifiersLocations,'stable');% I think stable is unneccessary because orderedIdenifiersLocations is sorted to begin with, but whatever...
end
orderOfIdentifiers2Use(length(orderedIdenifiersLocations)+1:end) = [];
recurrenceCounter = createRecurrenceCounter(orderedIdenifiersLocations);
maxIdentifierRecurrence = max(recurrenceCounter);
indexInMatrix = sub2ind([maxIdentifierRecurrence,numOfUniqueIdentifiers],recurrenceCounter,orderedIdenifiersLocations);
% orderOfIdentifiers2Use(isnan(indexInMatrix)) = [];    indexInMatrix(isnan(indexInMatrix)) = [];
end

function recurrenceCounter  = createRecurrenceCounter(vector)
% Gets a sorted vector and returns a result of the same size, where each
% location counts the number of (consecutive) recurrences of this value in vector.
efficientCount = true;
if efficientCount
    valueChangeIndecis = 1+find(diff(vector));
    recurrenceCounter = ones(size(vector));
    recurrenceCounter(valueChangeIndecis) = -(valueChangeIndecis-[1;valueChangeIndecis(1:end-1)]-1);
    recurrenceCounter = cumsum(recurrenceCounter);
else
    recurrenceCounter = ones(size(vector));
    for i=2:length(vector)
        if vector(i)==vector(i-1)
            recurrenceCounter(i) = recurrenceCounter(i-1)+1;
        end
    end
end
end