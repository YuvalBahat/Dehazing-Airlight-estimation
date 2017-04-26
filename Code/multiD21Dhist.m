function [histOutput1,histOutput2] = multiD21Dhist(x,histSecondArgument,normalizeOrExistingFigureHandle)
% Applying the matlab hist function on multidimensional data as if it was
% 1-dimensional (prevent the need to column stack). Besdies, setting the
% default bins number to 1000.
% Inputs:
% x - The (possibly) multi-dimensional data
% histSecondArgument - Any legitimate second argument to the matlab hist.
% normalize - 'norm' Display/Return a comulative normalized histogram. Pass
% [] in histSecondArgument for using default and returning normalized
% cumulative histogram  function. Passing 0 in histSecondArgument is equal
% to passing [] in histSecondArgument and 'norm' in normalize
% 
% Outputs:
% The same as in matlab's hist.
existingFigureHandle = [];
defaultBinsNum = false;
if nargin<2 || isempty(histSecondArgument)
    defaultBinsNum = true;
    if nargin<2
        normalizeOrExistingFigureHandle = 'yok';
    end
elseif histSecondArgument==0
    if nargin>2
        existingFigureHandle = normalizeOrExistingFigureHandle;
    end
    defaultBinsNum = true;
    normalizeOrExistingFigureHandle = 'norm';
elseif nargin<3
    normalizeOrExistingFigureHandle = 'yok';
end
if defaultBinsNum
    if any(isinf(x(:)) | imag(x(:))~=0);    error('Cannot automatically assign bins because of Inf or non-real values');    end;
    if any(round(x(:))-x(:))
        histSecondArgument = 1e3;
    else
        histSecondArgument = min(1e3,max(x(:))-min(x(:))+1);
    end
end

if iscell(x)
    if nargout~=0;  error('Not supporting cell array input with outputs');  end;
    displayCumulativeHist = true;
    allValsInX = [];    for i=1:length(x);  allValsInX = [allValsInX;x{i}(:)];  end; %#ok<AGROW>
    [~,histOutput2] = hist(allValsInX,histSecondArgument);
    histOutput1 = nan(length(x),histSecondArgument);
    for i=1:length(x)
        histOutput1(i,:) = hist(reshape(x{i},[],1),histOutput2);
    end
else
    [histOutput1,histOutput2] = hist(reshape(x,[],1),histSecondArgument);
    displayCumulativeHist = strcmp(normalizeOrExistingFigureHandle,'norm');
end
if displayCumulativeHist
    histOutput1 = bsxfun(@rdivide,cumsum(histOutput1.'),sum(histOutput1.'));
end

if nargout==0
    if displayCumulativeHist
        if ~isempty(existingFigureHandle)
            figure(existingFigureHandle);
            hold on;
            colorIndex = ceil(64*rand);
        else
            colorIndex = 1;
        end
        colors = colormap;
        plot(histOutput2,histOutput1,'Color',colors(colorIndex,:));
        if ~isempty(existingFigureHandle)
            hold off
        end
        ylabel('Portion of values');
    else
        bar(histOutput2,histOutput1);
    end
    clear histOutput1 histOutput2
end