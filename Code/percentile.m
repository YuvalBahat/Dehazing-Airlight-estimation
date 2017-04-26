function result = percentile(data,desiredPercentile)
perColumnPercentile = true;
if perColumnPercentile
    result = sort(data);
    if length(desiredPercentile)==3%Used in matrixFactorizationExperiment
        result = max(cat(1,desiredPercentile(2)*result(sub2ind(size(result),max(1,round(size(data,1)*desiredPercentile(1)))*ones(1,size(data,2)),...
            1:size(data,2))),result(sub2ind(size(result),max(1,round(size(data,1)*desiredPercentile(3)))*ones(1,size(data,2)),1:size(data,2)))),[],1);
    else
        if any(isnan(data(:)))
            perColumnIndex = round(sum(~isnan(data))*desiredPercentile);
        else
            perColumnIndex = round(size(data,1)*desiredPercentile)*ones(1,size(data,2));
        end
        result = result(sub2ind(size(result),max(1,perColumnIndex),1:size(data,2)));
    end
else
    result = sort(data(:));
    result = result(max(1,round(numel(data)*desiredPercentile)));
end