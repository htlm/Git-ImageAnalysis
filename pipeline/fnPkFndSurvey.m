%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Function to run fnPkFnd with all configurations of interest
%%%%%
%%%%% Last updated: 08/09/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pksArray = fnPkFndSurvey(expPath, overwriteImportedParams)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Determine peak parameters to import (later)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    overwriteImportedParams = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parse the path, determine filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathParse = strsplit(expPath,'/');
expName = char(pathParse(length(pathParse)));
post2017DatasetFlag = 0;

if isempty(expName)
    expName = char(pathParse(length(pathParse) - 1));
end

expIdentifier = char({expName(1:8)});

if str2double(expIdentifier(end-3:end)) > 2017
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
%%%%% Define Parameters List for peak finding and plate limits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if overwriteImportedParams
    parametersList = ...
        [120,   2800,   10;
        120,   3000,   9;
        150,   3000,   9;
        150,   5000,   9;
        160,   3000,   10;
        170,   800,    10;
        170,   2000,   10;
        170,   4000,   10;
        170,   5000,   10];
else
    parametersList = csvread('peakParams.csv', 0,0);
end

plateParams.radius = 525; 
plateCenters       = csvread('canaryPlateCenters.csv', 1, 0);
numImagesMatrix          = csvread('numImages.csv', 0, 0); 

%%%%% Assumes the input-type from plateCenters and numImages is the same as 
%%%%% str2num output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(str2double(expIdentifier) == plateCenters(:,1))
    plateParams.center = [525 , 555];
else
    I = str2double(expIdentifier) == plateCenters(:,1);
    
    plateParams.center = plateCenters(I,2:3);
end

if max(str2double(expIdentifier) == numImagesMatrix(:,1)) == 0
    plateParams.numImages = 288;
else
    I = str2double(expIdentifier) == numImagesMatrix(:,1);
    
    plateParams.numImages = numImagesMatrix(I,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Define Parameters List for peak finding and plate limits, assuming
%%%%% there are enough images to access the frames required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pksArray = {};

if plateParams.numImages > max(parametersList(:,1))
    for i = 1:size(parametersList,1)
        
        paramArray = parametersList(i,:);
        
        pksArray{i} = fnRunPkFnd(filenameArray, expPath, paramArray, plateParams);
        
    end
else
    'Not enough image frames'
end

end