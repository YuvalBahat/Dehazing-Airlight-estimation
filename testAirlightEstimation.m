% Estimating the airlight for an image in the 'images' sub-folder. Then
% displaying the estimated airlight, along with the image. 
% 
% A ground truth airlight can also be displayed for comparison. The GT airlight can be
% manually extracted by the user, or it can be taken from the corresponding
% file in the GTairlights sub-folder, if such file exists.

imageName = 'cones.png';
reextractAirlight = false;   % If true, manually extract the airlight color. Otherwise, use the pre-determined color, if exists.

addpath(genpath(pwd));
image2dehaze = im2double(imread(imageName));
estimatedAirlight = AirlightUsingPatchRecurrence(image2dehaze);

%% Displaying the estimated Airlight (and Ground truth airlight, if exists)

if reextractAirlight
    [GTairlight,GTairlightCovariance] = manuallyExtractAirlight(image2dehaze);
    GTairlightExists = true;
else
    GTairlightFileName = strcat('GTairlight_',imageName(1:end-4),'.mat');
    GTairlightExists = exist(GTairlightFileName,'file');
    load(GTairlightFileName)
end
figure;
subplot(4+GTairlightExists,1,1); imshow(createImageOfColor(estimatedAirlight,round(partialSize(image2dehaze,1:2).*[0.2,1])));
title(sprintf('Estimated Airlight:\n[%.2f,%.2f,%.2f]',estimatedAirlight));
subplot(4+GTairlightExists,1,2:4);   imshow(image2dehaze);   title('Hazy image')
if GTairlightExists
    subplot(5,1,5); imshow(createImageOfColor(GTairlight,round(partialSize(image2dehaze,1:2).*[0.2,1])));
    title(sprintf('GT Airlight:\n[%.2f,%.2f,%.2f]',GTairlight));
end