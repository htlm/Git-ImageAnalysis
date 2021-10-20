%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Create rgb images from background corrected.mat files. Assumes
%%%%%% saving of the new images in folders for the red-corrected-tiff,
%%%%%% green-corrected-tiff, rgb-correected-tiff (uses the entire expName
%%%%%% in the filename) and rgb-corrected-movie that has simpler filenames
%%%%%% for ffmpeg to use with a single wildcard. 
%%%%%%
%%%%%% Last updated: 08/09/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fnCanaryRGB(expPath, contrastVals, imagesToSkip)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the path, determine filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathParse = strsplit(expPath,'/');
expName = char(pathParse(length(pathParse)));
post2017DatasetFlag = 0;

if isempty(expName)
    expName = char(pathParse(length(pathParse) - 1));
end

expIdentifiers = char({expName(1:8)});

if str2double(expIdentifiers(end-3:end)) > 2017
    post2017DatasetFlag = 1;
end

redPath   = strcat(expPath, '/red-corrected/');
greenPath = strcat(expPath, '/green-corrected/');

if post2017DatasetFlag
    filenameArray = {strcat(redPath, expName, '_exc_DsRed_em_DsRed_channel2_t0');...
    strcat(greenPath, expName, '_exc_GFP_em_GFP_channel1_t0')};
else
    filenameArray = {strcat(redPath, expName, '_exc_DsRed_em_DsRed_t0');...
    strcat(greenPath, expName, '_exc_GFP_em_GFP_t0')};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image series basics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    imagesToSkip = [];
    if nargin < 2
        contrastVals = [];
    else
        maxRed      = contrastVals(1);
        maxGreen    = contrastVals(2);
    end
end

numImagesMatrix          = csvread('numImages.csv', 0, 0);

if max(str2double(expIdentifiers) == numImagesMatrix(:,1))==0
    numImages = 288;
    imageRange = [1 288];
else
    I = str2double(expIdentifiers) == numImagesMatrix(:,1);
    numImages = numImagesMatrix(I,2);
    imageRange = [1 numImages+length(imagesToSkip)];
end

startIm = imageRange(1);

if startIm < 10
        startImString = strcat('00', num2str(startIm));
else; if startIm < 100
            startImString = strcat('0', num2str(startIm));
      else
            startImString = num2str(startIm);
      end
end

data    = strcat(filenameArray{1}, startImString, '.mat');
load(data);
imageDimensions = size(hRed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocate matricies for .mat file data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp   = uint16(zeros(imageDimensions(1), imageDimensions(2), numImages));
expG  = uint16(zeros(imageDimensions(1), imageDimensions(2), numImages));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load matfiles - doesn't assume files are continuous
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = imageRange(1):imageRange(2)
    if isempty(imagesToSkip) || isempty(imagesToSkip == i) % If this image is NOT being skipped
        if i < 10
            dataName  = strcat(filenameArray{1}, '00', num2str(i), '.mat');
            dataNameG = strcat(filenameArray{2}, '00', num2str(i), '.mat');
        else; if i < 100
                dataName  = strcat(filenameArray{1}, '0', num2str(i), '.mat');
                dataNameG = strcat(filenameArray{2}, '0', num2str(i), '.mat');
              else
                dataName  = strcat(filenameArray{1}, num2str(i), '.mat');
                dataNameG = strcat(filenameArray{2}, num2str(i), '.mat');
              end
        end
        load(dataName);
        load(dataNameG);
        exp(:,:,i-startIm+1) = uint16(hRed);
        expG(:,:,i-startIm+1) = uint16(hGreen);
        clear dataName dataNameG hRed hGreen
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine contrast values, rescale and save as .tiff files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxRed   = double(max(max(max(exp))));
maxGreen = double(max(max(max(expG))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rescale images and save as .tiff files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expDouble  = double(exp);
expDoubleG = double(expG);

expRescale  = uint16(expDouble  ./ maxRed   .* 2^16); 
expRescaleG = uint16(expDoubleG ./ maxGreen .* 2^16); 

clear exp expG expDouble expDoubleG

redTiffPath   = strcat(expPath, '/red-corrected-tiff/');
greenTiffPath = strcat(expPath, '/green-corrected-tiff/');
rgbTiffPath = strcat(expPath, '/rgb-corrected-tiff/');
rgbTiffPathM = strcat(expPath, '/rgb-corrected-movie/');

if post2017DatasetFlag
    filenameArrayTiff = {strcat(redTiffPath, expName, '_exc_DsRed_em_DsRed_channel2_t0');...
    strcat(greenTiffPath, expName, '_exc_GFP_em_GFP_channel1_t0'); ...
    strcat(rgbTiffPath,   expName, '_rgb_t0'); ...
    strcat(rgbTiffPathM, 'rgb_t0')};
else
    filenameArrayTiff = {strcat(redTiffPath, expName, '_exc_DsRed_em_DsRed_t0');...
    strcat(greenTiffPath, expName, '_exc_GFP_em_GFP_t0'); ...
    strcat(rgbTiffPath,   expName, '_rgb_t0'); ...
    strcat(rgbTiffPathM, 'rgb_t0')};
end

% Setting saveastiff options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.overwrite   = 'true';
options.compression = 'no';
options.message     = false;

optionsRGB = options;
optionsRGB.color = true;

% Save files (red tiff, green tiff, rgb tiff with color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = imageRange(1):imageRange(2)
    if isempty(imagesToSkip) || isempty(imagesToSkip == i) % If this image is NOT being skipped
        if i < 10
            dataNameT     = strcat(filenameArrayTiff{1}, '00', num2str(i), '.tif');
            dataNameTG    = strcat(filenameArrayTiff{2}, '00', num2str(i), '.tif');
            dataNameTrgb  = strcat(filenameArrayTiff{3}, '00', num2str(i), '.tif');
            dataNameTrgbM = strcat(filenameArrayTiff{4}, '00', num2str(i), '.tif');
        else; if i < 100
                dataNameT     = strcat(filenameArrayTiff{1}, '0', num2str(i), '.tif');
                dataNameTG    = strcat(filenameArrayTiff{2}, '0', num2str(i), '.tif');
                dataNameTrgb  = strcat(filenameArrayTiff{3}, '0', num2str(i), '.tif');
                dataNameTrgbM = strcat(filenameArrayTiff{4}, '0', num2str(i), '.tif');
              else
                dataNameT     = strcat(filenameArrayTiff{1}, num2str(i), '.tif');
                dataNameTG    = strcat(filenameArrayTiff{2}, num2str(i), '.tif');
                dataNameTrgb  = strcat(filenameArrayTiff{3}, num2str(i), '.tif');
                dataNameTrgbM = strcat(filenameArrayTiff{4}, num2str(i), '.tif');
              end
        end
        hRed    = expRescale(:,:,i-startIm+1);
        hGreen  = expRescaleG(:,:,i-startIm+1);
        hRGB(:,:,1) = hRed;
        hRGB(:,:,2) = hGreen;
        hRGB(:,:,3) = uint16(zeros(size(hRed,1), size(hRed,2)));
        saveastiff(hRed, char(dataNameT), options);
        saveastiff(hGreen, char(dataNameTG), options);
        saveastiff(hRGB, char(dataNameTrgb), optionsRGB);
        saveastiff(hRGB, char(dataNameTrgbM), optionsRGB);
        clear hRed hGreen hRGB dataNameT dataNameTG dataNameTrgb
    end
end
end