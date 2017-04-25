imageFiles = dir(fullfile('images','*.png'));
imageFiles = cat(1,imageFiles,dir(fullfile('images','*.bmp')));
imageFiles = cat(1,imageFiles,dir(fullfile('images','*.jpg')));
clear estimatedAirlights
for fileNum = 1:length(imageFiles)
    startTime = now;
    estimatedAirlights(fileNum,:) = AirlightUsingPatchRecurrence(im2double(imread(fullfile('images',imageFiles(fileNum).name))));
    tocAndPrintWhenLarge(2,imageFiles(fileNum).name,(now-startTime)*24*3600);
end