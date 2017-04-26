function croppedImage = patchLikeCropAnImage(orgImage,patchSize)
% Cropping an image so that each pixel in croppedImage is the center of the patch extracted
% from orgImage.
if numel(patchSize)==1
    patchSize = patchSize*ones(1,2);
end
croppedImage = orgImage(floor(patchSize(1)/2)+1:end-floor(patchSize(1)/2),floor(patchSize(2)/2)+1:end-floor(patchSize(2)/2),:,:);