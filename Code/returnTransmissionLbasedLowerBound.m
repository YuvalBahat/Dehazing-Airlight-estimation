function [lowerBound,minimizingChannel] = returnTransmissionLbasedLowerBound(columnedImage,Airlight,Lle1ConstrainT,omega)
maxAerror = 0.1;
if ~any(size(Airlight,1)==[1,size(columnedImage,1),size(columnedImage,4)])
    error('Unsupported number of airlights');
end
if ~exist('omega','var')
    omega = 1;
end
if ~exist('Lle1ConstrainT','var')
    Lle1ConstrainT = false;
end
% Airlight = Airlight(1,:);
if ndims(columnedImage)>4
    error('Unsupported');
end
if ndims(Airlight)>2
    error('Didn''t check that yet since eliminating channels in which A is too low');
end
% if ~ismatrix(RGBimage)
%     columnedImage = reshape(RGBimage,[],size(RGBimage,3));
% else
%     columnedImage = RGBimage;
% end
% LgreaterThanZeroBound = max(1-bsxfun(@rdivide,max(eps,omega*columnedImage),max(eps,round(255*Airlight)/255)),[],2);
% [lowerBound,minimizingChannel] = max(1-bsxfun(@rdivide,max(eps,omega*columnedImage),max(eps,round(255*Airlight)/255)),[],2);
% Discarding channels in which A is too small for the lower bound to be
% reliable:
invalidChannels = Airlight<0.1;
validLessAirlights = find(all(invalidChannels,2));
if ~isempty(validLessAirlights)
    [~,chosenChannel] = max(invalidChannels(validLessAirlights,:),[],2);
    invalidChannels(sub2ind(size(invalidChannels),validLessAirlights,chosenChannel)) = false;
end
% [lowerBound,minimizingChannel] = max(1-bsxfun(@rdivide,max(eps,omega*columnedImage(repmat(validChannels,[1,1,partialSize(columnedImage,3:4)]))),...
%     Airlight(validChannels)),[],2);
oneMinusIoverA = 1-bsxfun(@rdivide,max(eps,omega*columnedImage),Airlight);
oneMinusIoverA(repmat(invalidChannels,[size(oneMinusIoverA,1)/size(Airlight,1),1,partialSize(oneMinusIoverA,3:4)])) = 0;
[lowerBound,minimizingChannel] = max(oneMinusIoverA,[],2);
% Now we assume A = A_true+Delta, where -maxAerror<=Delta<=maxAerror
if Lle1ConstrainT
    LsmallerThanOneBound = max(...
        min(cat(5,bsxfun(@rdivide,bsxfun(@minus,omega*columnedImage-maxAerror,Airlight),1-Airlight),...%For Delta <=0
        bsxfun(@rdivide,bsxfun(@minus,omega*columnedImage,Airlight),1-Airlight+maxAerror)),[],5),...%For Delta>=1-A_true
        [],2);%Adding a small additive term in the denominator to account for erroneous A
    LsmallerThanOneBound(any(columnedImage==1,2)) = 0; %Not using this bound when pixel is saturated, since L can be greater than 1 there.
    lowerBound = max(lowerBound,LsmallerThanOneBound);
end
lowerBound = max(lowerBound,0);
lowerBound(any(isnan(columnedImage),2)) = nan;
