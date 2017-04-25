function [executionInfo,varargout] = outOfMemoryResistantExcecution(executionStruct,varargin)
% executionStruct.functionName
% executionStruct.largeArgumentIndicator
% executionStruct.splittingInputDimension
% executionStruct.splittedOutputsIndicator
% executionStruct.splittingOutputDimension

functionName = executionStruct.functionName;
largeArgumentIndicator = executionStruct.largeArgumentIndicator;
splittingInputDimension = executionStruct.splittingInputDimension;
splittedOutputsIndicator = executionStruct.splittedOutputsIndicator;
splittingOutputDimension = executionStruct.splittingOutputDimension;

if isfield(executionStruct,'maxNumOfRetries')
    maxNumOfRetries = executionStruct.maxNumOfRetries;
else
    maxNumOfRetries = 4;
end
if isfield(executionStruct,'splittingFactor')
    splittingFactor = executionStruct.splittingFactor;
    maxNumOfRetries = min(maxNumOfRetries,length(splittingFactor));
else
    splittingFactor = 10.^(0:maxNumOfRetries-1);
end
problematicSize = size(varargin{find(largeArgumentIndicator,1,'first')},splittingInputDimension);
problematicSize = max(1,problematicSize); %For the trivial case of passing an empty argument.
if isfield(executionStruct,'trialNum2Start')
    trialNumber = executionStruct.trialNum2Start+1;
    previousBatchSize = floor(problematicSize/(splittingFactor(executionStruct.trialNum2Start)));
else
    trialNumber = 1;
end
% trialNumber = 1;    
success = false;
while ~success && maxNumOfRetries>trialNumber-1
    curBatchSize = floor(problematicSize/(splittingFactor(trialNumber)));
    if trialNumber>1
        if curBatchSize==previousBatchSize
            display(sprintf('Giving up after %d retrials due insufficient memory',trialNumber-2));
            rethrow(ME); %#ok<NODEF>
        end
    end
    curNumOfParts = ceil(problematicSize/curBatchSize);
    for partNum = 1:curNumOfParts
        executionString = [functionName,'('];
        for argInNum = 1:length(varargin)
            if curNumOfParts==1
                if ischar(varargin{argInNum})
                    executionString = [executionString,varargin{argInNum},','];
                else
                    executionString = [executionString,'varargin{',num2str(argInNum),'},'];
                end
            else
                if largeArgumentIndicator(argInNum)
                    argAssignmentstring = ['arg',num2str(argInNum),'=varargin{',num2str(argInNum),'}(',...
                        createArgumentIndexingSuffix(varargin{argInNum},splittingInputDimension)];
                    splittingDim = splittingInputDimension;
                    eval(argAssignmentstring);
                else
                    eval(['arg',num2str(argInNum),'=varargin{',num2str(argInNum),'};']);
                end
                if ischar(varargin{argInNum})
                    executionString = [executionString,varargin{argInNum},','];
                else
                    executionString = [executionString,'arg',num2str(argInNum),','];
                end
            end
        end
        executionString = [executionString(1:end-1),');'];
        outputsString = '[';
        for argOutNum = 1:nargout-1
            if splittedOutputsIndicator(argOutNum)
                outputsString = [outputsString,'argOut',num2str(argOutNum),','];
                if partNum==1
                    eval(['varargout{',num2str(argOutNum),'}=[];']);
                end
            else
                outputsString = [outputsString,'varargout{',num2str(argOutNum),'},'];
            end
        end
        executionString = [outputsString(1:end-1),'] = ',executionString];
        try
            eval(executionString);
            success = true;
            for argOutNum = 1:nargout-1
                if splittedOutputsIndicator(argOutNum)
                    eval(['varargout{',num2str(argOutNum),'}=cat(',num2str(splittingOutputDimension),...
                        ',varargout{',num2str(argOutNum),'},argOut',num2str(argOutNum),');']);
                end
            end
        catch ME
            switch ME.identifier
                case {'MATLAB:array:SizeLimitExceeded','MATLAB:nomem'}
                    success = false;
                    if maxNumOfRetries>trialNumber
                        display(sprintf('Out of memory in %s, trial #%d. Retrying...',functionName,trialNumber));
                    end
                    trialNumber = trialNumber+1;    previousBatchSize = curBatchSize;
                    break;
                otherwise
                    display(sprintf('outOfMemoryResistantExcecution does not support identifier %s',ME.identifier));
                    rethrow(ME);
            end
        end
    end
end
if ~success
    display(sprintf('Giving up after %d retrials due to insufficient memory',trialNumber-2));
    rethrow(ME); %#ok<NODEF>
end
executionInfo.numOfRetires = trialNumber-1;
executionInfo.finalBatchSize = curBatchSize;
executionInfo.numOfParts = curNumOfParts;
end

function indexingStringSuffix = createArgumentIndexingSuffix(variable,splittingDim)
indexingStringSuffix = [];
for dimNum = 1:ndims(variable)
    if dimNum==splittingDim
        indexingStringSuffix =...
            [indexingStringSuffix,'(partNum-1)*curBatchSize+1:',...
            'min(partNum*curBatchSize,size(varargin{argInNum},splittingDim)),'];
    else
        indexingStringSuffix = [indexingStringSuffix,':,'];
    end
end
indexingStringSuffix = [indexingStringSuffix(1:end-1),');'];
end