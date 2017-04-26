function [color,colorCovariances] = manuallyExtractAirlight(image)
h = figure;imshow((image));maximizeFigure;
title('Mark an area to sample color');
airlightRect = impoly;
airLightPixels = reshape(image(repmat(airlightRect.createMask,[1,1,3])),[],3);
color = mean(airLightPixels);
colorCovariances = cov(airLightPixels);
close(h);