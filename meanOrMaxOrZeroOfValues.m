function calculatedValues = meanOrMaxOrZeroOfValues(inputVals,dim,operation)
% Perform a min/max operation on an array, as long as the values on which
% we operate are not zero. If they are - returning zero. Used for cases
% where we ant to have the mean or the max value of an array, but cannot
% tolerate this array having zero values - in which case we would like to
% know that it has these values.

switch operation
    case 'mean'
        calculatedValues = mean(inputVals,dim).*~any(inputVals==0,dim);
    case 'max'
        calculatedValues = max(inputVals,[],dim).*~any(inputVals==0,dim);
    case 'min'
        calculatedValues = min(inputVals,[],dim).*~any(inputVals==0,dim);
    otherwise
        error('Operation %s unsupported',operation);
end