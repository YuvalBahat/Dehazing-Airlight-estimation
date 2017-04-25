function deisredSize = partialSize(A,deisredDimensions)
% Return the size of array A in the desired dimension deisredDimensions
Asize = size(A);
maxDesiredDim = max(deisredDimensions);
if maxDesiredDim>length(Asize)
    Asize(end+1:maxDesiredDim) = 1;
end
deisredSize = Asize(deisredDimensions);