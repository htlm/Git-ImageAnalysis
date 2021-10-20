%%%%% fnCanaryColonyProcessing - function to process all images for
%%%%% specified dataset into 3 Tidy Dataframe outputs (MATLAB tables) and
%%%%% two maps indicating pixel-colony allocation
%%%%%
%%%%%
%%%%% Last Updated: 05212020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [masterTable, colonyMask, Wmap, pkccMap, ccTable, metaData] ...
                = fnCanaryColonyProcessingNoT0(expPath, imageRange, saveFlag, optionalMaskFlag)
      tic      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Step 0: Set Parameters, load plate centers, import image
%%%%%         segementation information, parse filename info, 
%%%%%         peak finding parameters if applicable
%%%%% 
%%%%% Note:   This code CAN handle selective image selection if images are 
%%%%%         in an interval pattern. 
%%%%%
%%%%% Note:   This code CANNOT Handle selective frame removal. Images must
%%%%%         be contiguous or order must have a pattern!!!! 8/27/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathParse = strsplit(expPath,'/');
experiment = char(pathParse(length(pathParse)));
expIdentifier = experiment(1:8);
showMaskFlag = 0;

if str2double(expIdentifier(end-3:end)) > 2017
    post2017DatasetFlag = 1;
    if str2double(expIdentifier(end-3:end)) == 2019
        expFlag2019 = 1;
        maskThresh = 7000;
        
    else %%%  2018 dataset  %%%
        expFlag2019 = 0;
        maskThresh = 2500;
    end
else %%%  2017 dataset  %%%
    post2017DatasetFlag = 0;
    expFlag2019         = 0;
    maskThresh          = 2500;
end
colonyThresh = maskThresh*0.4;

%%%%% imageRange = [Start Image, Last Image, Image Interval, 
%%%%%               optional Image Interval] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metaData.numImages   = imageRange(2) - imageRange(1) + 1;
metaData.startIm     = imageRange(1);
metaData.imageRange  = imageRange;

plateParams.radius            = 525;
plateParams.erosion           = 1.1;%0.9;
plateParams.minPix            = 10;

if length(imageRange) > 2
    metaData.imageInterval = imageRange(3);
else
    metaData.imageInterval = 1;
end

%%%%% Access Plate Parameters and Centers #Tiago'sAnalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metaData.plateParams = plateParams;

plateCenters = csvread('canaryPlateCenters.csv', 1,0);

%%%%% Get Plate Center Info.    Assumes the input-type from 
%%%%% plateCenters.csv matches the str2double output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty( find(str2double(expIdentifier) == plateCenters(:,1)))

    plateParams.center = [525 , 555];
else
    I = find(str2double(expIdentifier) == plateCenters(:,1));
    if isempty(I)
        I = find(str2double(expIdentifier(2:end)) == plateCenters(:,1));
    end
    plateParams.center = plateCenters(I,2:3);
end

%%%%% Get Optimal Peak Finding Parameters.    Assumes the input-type from 
%%%%% Frame-Threshold-SZ.csv matches the str2double output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
peakParamDefaultList = csvread('peakParams.csv', 0,0);
peakParamsOptimized = csvread('Frame-Threshold-SZ.csv', 0,0);

I = find(str2double(expIdentifier) == peakParamsOptimized(:,1));
if peakParamsOptimized(I,2) == 0
    plateParams.peakParams = peakParamDefaultList;
else
    plateParams.peakParams = peakParamsOptimized(I,2:4);
end

%%%%% Load Peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if min(size(plateParams.peakParams)) == 1 %i.e. there is only one set of parameters/peaks
    pksName = strcat(expPath, '/pks_', num2str(plateParams.peakParams(1)), '-', ...
                                      num2str(plateParams.peakParams(2)), '-', ...
                                      num2str(plateParams.peakParams(3)), '.csv')
    pks = csvread(pksName, 0, 0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Step 1: Load images into 3D matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch char(expIdentifier)
    case '05232017'
        predata = strcat('CanaryData/', experiment, '/red-corrected/', ...
            '05232017_WT129_-7-6_NoQS_exc_DsRed_em_DsRed_t0');
        
        predataG = strcat('CanaryData/', experiment, '/green-corrected/', ...
            '05232017_WT129_-7-6_NoQS_exc_GFP_em_GFP_t0');
    otherwise
        if post2017DatasetFlag
            predata = strcat('CanaryData/', experiment, '/red-corrected/', ...
                experiment,'_exc_DsRed_em_DsRed_channel2_t0');
            
            predataG = strcat('CanaryData/', experiment, '/green-corrected/', ...
                experiment, '_exc_GFP_em_GFP_channel1_t0');
        else
            predata = strcat('CanaryData/', experiment, '/red-corrected/', ...
                experiment,'_exc_DsRed_em_DsRed_t0');
            
            predataG = strcat('CanaryData/', experiment, '/green-corrected/', ...
                experiment, '_exc_GFP_em_GFP_t0');
        end
end

data = strcat(predata, num2str(imageRange(2)), '.mat');
dataG = strcat(predataG, num2str(imageRange(2)), '.mat');

load(data, 'hRed');
load(dataG, 'hGreen');

hRedEnd   = uint16(hRed);
hGreenEnd = uint16(hGreen);

%%%%% Pre-allocate matricies for data storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp  = uint16(zeros(size(hRedEnd,1), size(hRedEnd, 2), ...
                    imageRange(2) - metaData.startIm + 1));
expG = uint16(zeros(size(hRedEnd,1), size(hRedEnd, 2), ...
                    imageRange(2) - metaData.startIm + 1));

intraColonyDist = zeros(size(hRedEnd,1), size(hRedEnd, 2), ...
                        imageRange(2) - metaData.startIm + 1);
spaceDeriv      = zeros(size(hRedEnd,1), size(hRedEnd, 2), ...
                        imageRange(2) - metaData.startIm + 1);
spaceDerivG     = zeros(size(hRedEnd,1), size(hRedEnd, 2), ...
                        imageRange(2) - metaData.startIm + 1);
timeDeriv       = zeros(size(hRedEnd,1), size(hRedEnd, 2), ...
                        imageRange(2) - metaData.startIm + 1);
timeDerivG      = zeros(size(hRedEnd,1), size(hRedEnd, 2), ...
                        imageRange(2) - metaData.startIm + 1);
dlogRdt         = zeros(size(hRedEnd,1), size(hRedEnd, 2), ...
                        imageRange(2) - metaData.startIm + 1);
dlogGdt         = zeros(size(hRedEnd,1), size(hRedEnd, 2), ...
                        imageRange(2) - metaData.startIm + 1);
                    
%%%%% Load Background Corrected Data (.mat files)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp(:,:,end) = hRedEnd;
expG(:,:,end) = hGreenEnd;

for i = metaData.startIm : metaData.imageInterval : metaData.imageRange(2)
    if i < 10
        dataName = strcat(predata, '00', num2str(i), '.mat');
        dataNameG = strcat(predataG, '00', num2str(i), '.mat');
    else if i < 100
            dataName = strcat(predata, '0', num2str(i), '.mat');
            dataNameG = strcat(predataG, '0', num2str(i), '.mat');
        else
            dataName = strcat(predata, num2str(i), '.mat');
            dataNameG = strcat(predataG, num2str(i), '.mat');
        end
    end
    load(dataName,  'hRed');
    load(dataNameG, 'hGreen');
    if i > 1
    exp(:,:,i-metaData.startIm+1)  = uint16(hRed) - exp(:,:,1);
    expG(:,:,i-metaData.startIm+1) = uint16(hGreen) - exp(:,:,1);
    else 
        continue % No data is put in the matrix for the first image in imageRange
    end
    clear dataName dataNameG hRed hGreen
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Step 2: Make (and Save) [or Load (if applicable)] the Mask. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 3
    basicMask = 1 - erodePlate(imbinarize(exp(:,:,end)), [526, 550], 48, 1, 0);
else
    basicMask = true(size(hRedEnd,1), size(hRedEnd, 2));
end

bwRedStack = false(size(hRedEnd,1), size(hRedEnd, 2), imageRange(2) - metaData.startIm + 1);

%%%%% Erode the plate edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
erodeMask = erodePlate(imbinarize(exp(:,:,end)), plateParams.center, ...
                       plateParams.radius, plateParams.erosion, 0);

%%%%% Create and save final mask - figure out these floating parameters!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:imageRange(2) - metaData.startIm + 1 
    if expFlag2019
        bwRedStack(:,:,i) = logical( bpass(double(exp(:,:,i)), 0, [], maskThresh))...
                        .* erodeMask.*basicMask;
    else
    bwRedStack(:,:,i) = logical( bpass(double(exp(:,:,i)), 0, [], maskThresh))...
                        .* erodeMask.*basicMask;
    end
end

bwRed = logical(sum(bwRedStack, 3));
colonyMask = bwRed;
metaData.numPix = sum(sum(colonyMask));
maskName = strcat('CanaryData/', experiment, '/', expIdentifier, '_endMask.mat');
save(maskName, 'bwRed');
expMask = repmat(bwRed, [1, 1, metaData.numImages]);

if showMaskFlag
    imshow(bwRed)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Step 3: Create and Save colony pixel-allocation maps for each set of 
%%%%%         peaks identified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Dilate the peak locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M   = zeros(size(hRedEnd));
M(sub2ind(size(hRedEnd),pks(:,2),pks(:,1))) = 1;
MM  = imdilate(M,ones(4));

%%%%% Impose minima in the peak locations on the colony mask 
%%%%% Using the mask as a baseline because it overwrites phage bursts
%%%%% Use the mask dilated with minima at the peak locations to calculate
%%%%% the watershed and get the final pixel-to-colony allocations in
%%%%% variable WMap. (Saved)
%%%%%
%%%%% Assumption: There are never more than 2^16-1 colonies/plate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = imimposemin((2^16 - 1) - hRedEnd, MM);
W = watershed(double(N).*double(bwRed));
preWmap = (uint16(W) - 1).*uint16(bwRed);
W2 = bwareaopen(preWmap, plateParams.minPix);
Wmap = preWmap.*uint16(W2);

if showMaskFlag
    figure
    subplot(1,3,1)
    imagesc(W)
    caxis([0 100])
    title('Watershed')
    
    subplot(1,3,2)
    imagesc(W2)
    title('removing small objects?')
    
    subplot(1,3,3)
    imagesc(Wmap)
    title('Final pixel-to-colony map')
end
    
WmapName = strcat('CanaryData/', experiment, '/', expIdentifier, '_Wmap.mat');
save(WmapName, 'Wmap');
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Step 4: Pixel-based filtering in time [no filtering in spacial dim]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Time dimension filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp = movmean(exp, 5, 3);
expG = movmean(expG, 5, 3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Step 5: Create meshgrid and apply mask to create pixel-based columns.
%%%%%         Create basicTable!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Create Meshgrid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xArray = 1:1:size(hRedEnd, 2);
yArray = 1:1:size(hRedEnd, 1);
tArray = imageRange(1):metaData.imageInterval:imageRange(2);

[xBlock, yBlock, tBlock] = meshgrid(xArray, yArray, tArray);
clear xArray yArray tArray
xBlock = uint16(xBlock);
yBlock = uint16(yBlock);
tBlock = uint16(tBlock);

%%%%% Calculating Temporal and Spatial Derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:imageRange(2) - metaData.startIm + 1
    
    %%%%% Take the spatial derivative
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spaceDeriv(:,:,i)      = 4*del2(double(exp(:,:,i)));
    spaceDerivG(:,:,i)     = 4*del2(double(expG(:,:,i)));
    intraColonyDist(:,:,i) = bwdist((exp(:,:,i)<1000), 'euclidean');
    
    %%%%% Take the time derivative
    %%%%% Time interval = imageInterval*(1/6th hour)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i == imageRange(2) - metaData.startIm + 1
        continue
    else
        timeDeriv(:,:,i)     = diff(double(exp(:,:,i:i+1)), [], 3)...
                                ./ (double(metaData.imageInterval)*1/6); 
        timeDerivG(:,:,i)    = diff(double(expG(:,:,i:i+1)), [], 3)...
                                ./ (double(metaData.imageInterval)*1/6);
    
        dlogRdt(:,:,i)  = diff(log(double(exp(:,:,i:i+1))), [], 3)...
                                ./ (double(metaData.imageInterval)*1/6); 
        dlogGdt(:,:,i)  = diff(log(double(expG(:,:,i:i+1))), [], 3)...
                                ./ (double(metaData.imageInterval)*1/6);
    end
end

%%%%% Basic independent variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x     = xBlock(expMask);
y     = yBlock(expMask);
t     = double(tBlock(expMask)).*10/60; 
tInt  = tBlock(expMask);
tCtgl = categorical(round(t, 3)); 
expID = categorical(cellstr(repmat(expIdentifier, length(x), 1))); 


%%%%% Pixel-based intensity/spreading information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = double(exp(expMask));
logR = log(double(exp(expMask)));
G = double(expG(expMask));
logG = log(double(expG(expMask)));
dTime = timeDeriv(expMask);
dTimeG = timeDerivG(expMask);
dlogRdt = timeDeriv(expMask);
dlogGdt = timeDerivG(expMask);
dSpace = spaceDeriv(expMask);
dSpaceG = spaceDerivG(expMask);

%%%%% Pixel-based location information - VERIFY LATER!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intraColonyDistR = intraColonyDist(expMask);
distToEdge  = double(plateParams.radius) - ...
              sqrt((double(plateParams.center(1)) - double(x)).^2 + ...
                   (double(plateParams.center(2)) - double(y)).^2);
invDistToEdge  = 1./distToEdge;
polarX = double(x) - double(plateParams.center(1));
polarY = double(y) - double(plateParams.center(2));
[theta, radius] = cart2pol(polarX, polarY);

%%%%% Create basicTable - smaller version of masterTable before colony info
%%%%% is added
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metaData.vars       = {'expID', 'x', 'y', 't', 'tInt', 'tCtgl', 'R', 'G', 'logR', 'logG', ...
                 'distToEdge', 'invDistToEdge', 'theta', 'radius', ...
                 'intraColonyDistR', 'dTime', 'dTimeG', 'dlogRdt', 'dlogGdt', ...
                 'dSpace', 'dSpaceG'};
metaData.varTypes   = {'categorical', 'uint16', 'uint16', 'double', 'uint16', ...
                       'categorical', 'double', 'double', 'double', 'double', ...
                 'double', 'double', 'double', 'double', ...
                 'double', 'double', 'double', 'double', 'double', ...
                 'double', 'double'};
             
basicTable = table(expID,  x,   y,   t,   tInt,   tCtgl,   R,   G,   logR,   logG,  ...
                  distToEdge,   invDistToEdge,   theta,   radius,  ...
                  intraColonyDistR,   dTime,   dTimeG,   dlogRdt,   dlogGdt,  ...
                  dSpace,   dSpaceG, ...
                  'VariableNames', metaData.vars); 
              
tableToUse = basicTable;

%%%% Save basicTable as a baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFlag
    save(strcat(expIdentifier, 'basicTable-NoT0-',  date, '.mat'), 'basicTable', '-v7.3');
else
    'Not saving basicTable'
end
'basicTableSaved'
clear basicTable

clear exp expG R logR G logG dTime dTimeG dlogRdt dlogGdt dSpace dSpaceG ...
      intraColonyDistR distToEdge invDistToEdge polarX polarY theta radius


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Step 6: Calculate colony-based data, create ccTable, 
%%%%%         Add all data together to create masterTable!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
%%%%% Columns based on Connected Component (Colony) information
%%%%%   I am going to *assume* that since the connected components and
%%%%%   regionprops info are coming from the same Wmap, that the order of
%%%%%   the components in the output is THE SAME. (2/25/18)
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Combining Steps 6 and 7 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Step 7: Add all data together to create basicTableW = basicTable with
%%%%%         the watershed information included/propogated from ccTable
%%%%%
%%%%% Variables to be transitioned from ccTable: 
%%%%%         colonyID, area, totalColonyIntR, totalColonyIntG
%%%%%
%%%%% Note: Not all pixels in basicTable may be in one of the final
%%%%%       connected component pixel lists. CHECK THIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

regProps = regionprops(Wmap, hRedEnd, ...
    'pixellist','PixelValues','MajorAxisLength', 'MinorAxisLength', ...
    'Centroid', 'Area', 'MeanIntensity');
regPropsTable = struct2table(regProps);
regPropsTable.ccID = [1:size(regProps,1)]';
regPropsTable = regPropsTable(regPropsTable.Area >= 50, :);

%%%%% Build an association between the connected components and the 
%%%%% colony peaks (using data from Wmap): pkccMap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pkccMap  = [pks ones(length(pks),1)];
pkccMap  = double(pkccMap);
mapCheck = zeros(length(pks(:,1)), 1);
for j = 1:length(pks)
    
    pkccMap(j,3) = Wmap(pks(j,2), pks(j,1)); 
    %%% Verify that the peak is contained in the cc identified by Wmap, or
    %%% 0 in Wmap
    ccIdx = find(regPropsTable.ccID == pkccMap(j,3));
    if ~isempty(ccIdx)
    mapCheck(j) = 0 < find( and( pks(j,1) == regPropsTable.PixelList{ccIdx, 1}(:,1), ...
                                 pks(j,2) == regPropsTable.PixelList{ccIdx, 1}(:,2)));
    else %%% A peak may have had its region removed if it was a small object
        mapCheck(j) = 2;
    end
end
if isempty(find(mapCheck == 2)) 
    'All small object peaks removed'
else
    'Small object peak found'
    smallPks                = find(mapCheck == 2);
    pks(smallPks, :)        = [];
    pkccMap(smallPks, :)    = [];
end
    

if ~isempty(find(mapCheck == 0))
    'Error: CCID MisMatch!!'
end

%%%%% Calculate the number of pixels to be included in the final table
%%%%% based on the connected component information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tableNumPix           = zeros(height(regPropsTable), 1);

for i = 1:height(regPropsTable)
    tableNumPix(i) = length(regPropsTable.PixelList{i});
end

tableNumPix = tableNumPix(tableNumPix ~= 0); % Shouldn't change anything, all colonies with 0 pix removed above
finalTableHeight = sum(tableNumPix) * metaData.numImages;
metaData.colonyNumPix = tableNumPix;

%%%%% Pre-allocate colony-centric data matricies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metaData.varsExt      = [metaData.vars, {'colonyID', 'colonyIDctgl', 'area', ...
                        'totalColonyIntR', 'totalColonyIntG'}];
metaData.varTypesExt  = [metaData.varTypes, {'uint16', 'categorical', 'uint16', ...
                        'double', 'double'}];
basicTableW           = table('Size', ...
                        [finalTableHeight, ...
                         length(metaData.varTypesExt)], ...
                        'VariableTypes', metaData.varTypesExt, ...
                        'VariableNames', metaData.varsExt);

colonyArea            = zeros(metaData.numImages,length(pks));
totalColonyIntensityR = zeros(metaData.numImages,length(pks));
totalColonyIntensityG = zeros(metaData.numImages,length(pks));

metaData.pksToRegPropsLinker   = zeros(height(regPropsTable), 1);

timeToHarvestPix      = zeros(height(regPropsTable), 1);
timeToInsertPix       = zeros(height(regPropsTable), 1);

    toc
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Step 8: Iterate over peaks to identify the corresponding connected component
%%%%% and create the subTable cooresponding to that cc
%%%%% 
%%%%% Future upgrade [Optional]: iterate over the connected components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(pks)

    pkOI    = pks(i,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%   Isolate the connected component containing the peak of interest
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ccOI      = pkccMap(and(pkccMap(:,1)==pkOI(1), pkccMap(:,2)==pkOI(2)),3);

    regPropsTableIdx = find(regPropsTable.ccID == ccOI);
    metaData.pksToRegPropsLinker(i) = regPropsTableIdx;
    ccPixList = regPropsTable.PixelList{regPropsTableIdx};
    
    %%%%% Preallocate R and G arrays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ccPixRArray = zeros(metaData.numImages, regPropsTable.Area(regPropsTableIdx));
    ccPixGArray = zeros(metaData.numImages, regPropsTable.Area(regPropsTableIdx));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Find the indicies of the basicTable that have the pixels of 
    %%%%% interest for each cc and builds colonyTable that includes all
    %%%%% pixel information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    colonyTable = table('Size', ...
                        [regPropsTable.Area(regPropsTableIdx) * metaData.numImages, ...
                         length(metaData.varTypes)], ...
                        'VariableTypes', metaData.varTypes, ...
                        'VariableNames', metaData.vars);
                    
    for j = 1:regPropsTable.Area(regPropsTableIdx) % Adds each pixel to the subset of masterTable
        pixOI = find( and( regPropsTable.PixelList{regPropsTableIdx}(j,1) == tableToUse.x,...
                           regPropsTable.PixelList{regPropsTableIdx}(j,2) == tableToUse.y)  );        
        
        colonyTable((j-1)*(metaData.numImages) + 1 : j*(metaData.numImages), :) = tableToUse(pixOI, :);
        
        ccPixRArray(:,j)  = table2array(tableToUse(pixOI, 'R')); 
        ccPixGArray(:,j)  = table2array(tableToUse(pixOI, 'G')); 
    end
    timeToHarvestPix(i) = toc;
    
    
    %%%%% Threshold the background to identify colony-specific info
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    backgroundThresh = colonyThresh; %%%%% ARBITRARY VALUE - is checked later in downstream processing
    tccPixArrayL     = ccPixRArray > backgroundThresh;
    
    %%%%% Calculate total colony Area and Intensity as f(t)
    %%%%% Take the sum row by row (each row is a single timepoint)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    colonyArea(:,i)              = sum(tccPixArrayL,2); 
    totalColonyIntensityR(:,i)   = sum(ccPixRArray.*tccPixArrayL,2);
    totalColonyIntensityG(:,i)   = sum(ccPixGArray.*tccPixArrayL,2); 

    
    %%%%% Duplicate the colony-specific data to the size of the
    %%%%% corresponding number of pixels for each colony in the colonyTable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    colonyTable                  = sortrows(colonyTable, {'t', 'x'});
    colonyTable.colonyID         = ccOI .* uint16(ones(height(colonyTable),1));
    colonyTable.colonyIDctgl     = categorical(colonyTable.colonyID);
    areaMatrix                   = repmat(uint16(colonyArea(:,i)), 1, length(ccPixList(:,1)))';
    colonyTable.area             = areaMatrix(:);
    totalColonyIntRMatrix        = repmat(totalColonyIntensityR(:,i), 1, length(ccPixList(:,1)))';
    colonyTable.totalColonyIntR  = totalColonyIntRMatrix(:);
    totalColonyIntGMatrix        = repmat(totalColonyIntensityG(:,i), 1, length(ccPixList(:,1)))';
    colonyTable.totalColonyIntG  = totalColonyIntGMatrix(:);

    basicTableW(sum(tableNumPix(1:regPropsTableIdx-1))*metaData.numImages +1 : ...
                sum(tableNumPix(1:regPropsTableIdx-1))*metaData.numImages + tableNumPix(regPropsTableIdx)*metaData.numImages, :) = colonyTable;

    timeToInsertPix(i) = toc;
end
toc
clear colonyTable colonyArea areaMatrix totalColonyIntensityR totalColonyIntRMatrix ...
      totalColonyIntensityG totalColonyIntGMatrix ccPixRArray ccPixGArray

metaData.timeToHarvestPix = timeToHarvestPix;
metaData.timeToInsertPix  = timeToInsertPix;

%%%%% Last table edits - categorical variable and addition of reference
%%%%% categories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove extra pre-allocated rows
basicTableW = basicTableW(basicTableW.t ~=0, :);

expIDmin = min(str2double(char(unique(basicTableW.expID))));
if expIDmin < 1e7
    expIDminstr = strcat('0', num2str(expIDmin));
else
    expIDminstr = num2str(expIDmin);
end

basicTableW.tCtgl = addcats(basicTableW.tCtgl, 'ref', 'before', ...
                            num2str(round(min(basicTableW.t), 3)));
basicTableW.colonyIDctgl = addcats(basicTableW.colonyIDctgl, 'ref', ...
                                   'before', char(num2str(min(basicTableW.colonyID))));
basicTableW.expID = addcats(basicTableW.expID, 'ref', 'before', expIDminstr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Step 9: Calculate local density information 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Calculate Local Density columns, add to basicTable
%%%%% Upgrade [Optional]: Calculate for each pixel, not just peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tableToUse = basicTableW;

localDensityThreshMatrix = 25:10:150;
[localDensityTable, metaData] = fnLocalDensityTableCreation(tableToUse, metaData, pks, pkccMap, ...
                                                localDensityThreshMatrix);
toc                                          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Step 10: Save masterTable, pixelTable, ccTable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tableToUse = localDensityTable;
masterTable = tableToUse; 

%%%%% Generating last metadata parts
%%%%% Upgrade [Optional]: Add data on muMax, spatial distribution, etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metaData.maxRed         = max(masterTable.R);
metaData.maxGreen       = max(masterTable.G);
metaData.numColonies    = length(unique(masterTable.colonyID));
metaData.hRedEnd        = hRedEnd;

%%%%% Create ccTable here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ccTable               = table('Size', ...
                        [size(pks, 1) * metaData.numImages, ...
                         length(metaData.varTypesExt2)], ...
                        'VariableTypes', metaData.varTypesExt2, ...
                        'VariableNames', metaData.varsExt2);

timeToAddCC      = zeros(height(regPropsTable), 1);

for i = 1:length(pks)
    pkOI    = pks(i,:);
    pkTable = tableToUse(and(tableToUse.x==pkOI(1), tableToUse.y==pkOI(2)),:); 
        
    ccTable((i-1)*(metaData.numImages) + 1 : i*(metaData.numImages), :) ...
            = pkTable;

    timeToAddCC(i) = toc;
end

toc
metaData.timeToAddCC = timeToAddCC;

%%%%% Saving Tables and Metadata. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFlag
    'Saving...'
    %%%% basicTable was saved and cleared above
    save(strcat(expIdentifier, 'masterTable-NoT0-', date, '.mat'), 'masterTable', '-v7.3');
    save(strcat(expIdentifier, 'ccTable-NoT0-',     date, '.mat'), 'ccTable', '-v7.3');
    save(strcat(expIdentifier, 'pkccMap-NoT0-',     date, '.mat'), 'pkccMap', '-v7.3');
    save(strcat(expIdentifier, 'metaData-NoT0-',    date, '.mat'), 'metaData', '-v7.3');
    
    metaDataName = strcat('CanaryData/', experiment, '/', expIdentifier, ...
                          '_metaData.mat');
    save(metaDataName, 'metaData');
else
    'Not Saving'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Calculates local density information for a set of local density
    %%%%% thresholds
    %%%%%
    %%%%% Last Updated: 02/18/2019
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [localDensityTable, metaData] ...
            = fnLocalDensityTableCreation(tableToUse, metaData, pks, pkccMap, localDensityThreshMatrix)
        
        %%%%% Plan the variables types and names to be calculated
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        newVarTypes = cell(1, length(localDensityThreshMatrix) + 1);
        newVarTypes(:) = {'uint16'};
        metaData.varTypesExt2 = [metaData.varTypesExt, newVarTypes];
        
        newVarNames = cell(1, length(localDensityThreshMatrix) + 1);
        newVarNames(1) = {'hundredPix'};
        for k = 1:length(localDensityThreshMatrix)
            localDensityThresh = localDensityThreshMatrix(k);
            
            %%%%% Generate and Store LDName
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            LDName = strcat('nr', num2str(localDensityThresh));
            newVarNames(k+1) = {LDName};
        end
        
        metaData.varsExt2 = [metaData.varsExt, newVarNames];
        
        %%%%% Pre-allocate the table
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        localDensityTable = table('Size', ...
            [height(tableToUse), ...
            length(metaData.varTypesExt2)], ...
            'VariableTypes', metaData.varTypesExt2, ...
            'VariableNames', metaData.varsExt2);
        
        %%%%% Generate Local Density Data (and store variable names in the metaData)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pixelToPks = squareform(pdist(pks));
        
        for k = 1:length(pks)
            pkOI    = pks(k,:); %%% Each peak has a connected component associated with it from earlier in the algorithm
                
                ccOI        = pkccMap(and(pkccMap(:,1)==pkOI(1), pkccMap(:,2)==pkOI(2)),3);
                metaDataIdx = metaData.pksToRegPropsLinker(k);
                
                if ccOI == 0
                    'Error: ccOI == 0'
                    continue
                else
                    subTable  = tableToUse(tableToUse.colonyID == ccOI,:);
                    subTable.hundredPix = repmat(sum(pixelToPks(k,:) < 100), ...
                        length(subTable.t), 1);
                    for l = 1:length(localDensityThreshMatrix)
                        
                        localDensityThresh = localDensityThreshMatrix(l);
                        LDName = newVarNames{l+1};
                        
                        subTable.(LDName) = repmat(sum(pixelToPks(k,:) < localDensityThresh), ...
                            length(subTable.t), 1);
                        clear LDName
                    end
                    
                    localDensityTable(sum(metaData.colonyNumPix(1:metaDataIdx-1))*metaData.numImages +1 : ...
                        sum(metaData.colonyNumPix(1:metaDataIdx-1))*metaData.numImages ...
                        + metaData.colonyNumPix(metaDataIdx)*metaData.numImages, :) ...
                        = subTable;
                end
        end
    end


end