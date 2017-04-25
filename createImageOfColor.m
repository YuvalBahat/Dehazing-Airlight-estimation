function colorIm = createImageOfColor(color,imSize,frameWidth)
if nargin<3
    frameWidth = 0;
    if nargin<2
        imSize = 200;
    end
end
colorIm = bsxfun(@times,ones(imSize),reshape(color,[1,1,3]));
if frameWidth>0
    colorIm(1:frameWidth,:,:) = 0;
    colorIm(end-frameWidth+1:end,:,:) = 0;
    colorIm(:,1:frameWidth,:) = 0;
    colorIm(:,end-frameWidth+1:end,:) = 0;
end
if nargout<1
    figure;
    imshow(colorIm);
    clear colorIm
end