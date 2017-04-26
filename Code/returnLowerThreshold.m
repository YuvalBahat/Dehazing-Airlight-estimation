function lowerThreshold = returnLowerThreshold(values,desiredNumRange,fixedMinValue)
% Returning a threshold value satisfying the following constraints (when desiredNumRange>1, meaning not dealing with portion of values):
% 1.Minimum of desiredNumRange(1) passing the threshold
% 2.Above this minimum, all are above fixedMinValue
% 3.Maximum of desiredNumRange(2), if exists
lowerThreshold = sort(values,'descend');
if desiredNumRange(1)<=1 && desiredNumRange(1)~=0
    lowerThreshold = lowerThreshold(round(desiredNumRange(1)*length(lowerThreshold)));
else
    if length(desiredNumRange)==1; desiredNumRange(2) = length(values); end;
    if desiredNumRange(1)==0
        lowerThreshold = max(fixedMinValue,lowerThreshold(min(length(lowerThreshold),desiredNumRange(2))));
    else
        lowerThreshold = max(min(fixedMinValue,...
            lowerThreshold(min(length(lowerThreshold),desiredNumRange(1)))),lowerThreshold(min(length(lowerThreshold),desiredNumRange(2))));
    end
end