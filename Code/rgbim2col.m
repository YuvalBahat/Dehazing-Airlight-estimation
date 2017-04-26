function B = rgbim2col(A,patchSize,distinctOrSliding)
if strcmp(distinctOrSliding,'distinctWithoutPadding')
    A = A(1:floor(size(A,1)/patchSize(1))*patchSize(1),1:floor(size(A,2)/patchSize(2))*patchSize(2),:);
    distinctOrSliding = 'distinct';
end
B = im2col(A(:,:,1),patchSize,distinctOrSliding);
B = cat(3,B,nan([partialSize(B,1:2),size(A,3)-1]));
for channelNum = 2:size(A,3)
    B(:,:,channelNum) = im2col(A(:,:,channelNum),patchSize,distinctOrSliding);
end