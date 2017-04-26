% Estimating the airlight on all images in the 'images' sub-folder, and
% comparing to the GT airlight for those images that contain regions of
% sky.

imageFiles = dir(fullfile('images','*.png'));
imageFiles = cat(1,imageFiles,dir(fullfile('images','*.bmp')));
imageFiles = cat(1,imageFiles,dir(fullfile('images','*.jpg')));

estimatedAirlights = nan(length(imageFiles),3);
AirlightDistances = nan(length(imageFiles),1);

for fileNum = 1:length(imageFiles)
    startTime = now;
    GTfilename = ['GTairlight_',imageFiles(fileNum).name(1:end-4),'.mat'];
    estimatedAirlights(fileNum,:) = AirlightUsingPatchRecurrence(im2double(imread(fullfile('images',imageFiles(fileNum).name))));
    if exist(GTfilename,'file')
        load(GTfilename)
        AirlightDistances(fileNum) = sqrt(sum((estimatedAirlights(fileNum,:)-GTairlight).^2));
    end
    tocAndPrintWhenLarge(2,sprintf('Processed image %s',imageFiles(fileNum).name),(now-startTime)*24*3600);
end

display(sprintf('Average distance between estimated and GT airlights (in RGB space) is %.2f',nanmean(AirlightDistances)));