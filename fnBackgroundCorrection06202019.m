% Performing background correction on Canary images given the Typhoon
% Correction
%
% Last updated: 07/30/2018

% %% Create the transformation

function fnBackgroundCorrection06202019(filenameRed, filenameGreen, expPath, ...
                                       saveFlag, bkgndCorrYr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for saveFlag and path info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
    saveFlag = 1;

        if nargin < 3
            expPath = '';
        end
end

pathParse = strsplit(expPath,'/');
experiment = char(pathParse(length(pathParse)));
expIdentifier = experiment(1:8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set any arbitrary parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
canaryBkgnd1Red      = 1500; % Keeping the background subtraction used in 
                             % calculating the background correction 
canaryBkgnd1Green    = 1500; % Keeping the background subtraction used in 
                             % calculating the background correction 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load correction information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if post2017DatasetFlag
%     load('mdlRed.mat');
%     load('mdlGreen.mat');
% else
%     load('mdlRed2017-2.mat');
%     load('mdlGreen2017-2.mat');
% end

switch bkgndCorrYr
    
    case 2017
        load('mdlRed-2017_06202019.mat');
        load('mdlGreen-2017_06202019.mat');
        
    case 2018
        load('mdlRed-2018_06202019.mat');
        load('mdlGreen-2018_06202019.mat');
        
    case 2019
        load('mdlRed-2019_06202019.mat');
        load('mdlGreen-2019_06202019.mat');
        
end


coefsR = mdlRed.Coefficients.Estimate;
coefsG = mdlGreen.Coefficients.Estimate;

% Import images and subtract basic Canary background
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
canaryRedRaw   = imread(char(strcat(strcat(expPath, '/DsRed-DsRed/'), strcat(filenameRed, '.tif'))));
canaryGreenRaw = imread(char(strcat(strcat(expPath, '/GFP-GFP/'),   strcat(filenameGreen, '.tif'))));

canaryRedIm       = canaryRedRaw - canaryBkgnd1Red; 
canaryGreenIm      = canaryGreenRaw - canaryBkgnd1Green; 

canaryRed       = double(canaryRedIm);
canaryGreen     = double(canaryGreenIm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the colony information (Test which of the below work better)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
canaryMaskRed     = canaryRedIm   > 0;
% colonyMask = im2bw(imadjust(uint16(canaryRedIm)));
% redColoniesOrigCanary0 = double(canaryRedIm).*colonyMask;
% greenColoniesOrigCanary0 = double(canaryGreenIm).*colonyMask;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the correction and save output images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m, n]  = size(canaryRed);
iMatrixR = 1:m;
jMatrixR = 1:n;
[xR, yR] = meshgrid(iMatrixR, jMatrixR);
[xG, yG] = meshgrid(iMatrixR, jMatrixR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Distance To Center Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plateCenters = csvread('canaryPlateCenters2.csv', 0,0);

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
%%%% Looks good investigating this calculation in the 110062017 dataset
%%%% on 11302018
distToCenter = sqrt((xR - plateParams.center(2)).^2 + (yR - plateParams.center(1)).^2)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2017, 2018 and 2019 corrections now all have the same form: 
%  
%  Red Correction: typhoonColoniesRed ~ 1 + xDir + yDir + canaryColoniesRed ...
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   + distToCenter + ...
%                                      + xDir*yDir ...
%                                      + xDir*canaryColoniesRed ...
%                                      + yDir*canaryColoniesRed ...
%                                      + distToCenter*canaryColoniesRed;
xyCorrectionR =  coefsR(1) ...
    + coefsR(2).*xR' + coefsR(3).*yR' ...
    + coefsR(4).*canaryRed.*canaryMaskRed ...
    + coefsR(5).*distToCenter ...
    + coefsR(6).*xR'.*yR' ...
    + coefsR(7).*xR'.*canaryRed.*canaryMaskRed ...
    + coefsR(8).*yR'.*canaryRed.*canaryMaskRed ...
    + coefsR(9).*distToCenter.*canaryRed.*canaryMaskRed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2017, 2018 and 2019 corrections now all have the same form: 
%  
%  Green Correction: typhoonColoniesGreen ~ 1 + xDir + yDir + canaryColoniesGreen ...
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   + distToCenter ...
%                                          + xDir*yDir ...
%                                          + xDir*canaryColoniesGreen ...
%                                          + yDir*canaryColoniesGreen ...
%                                          + distToCenter*canaryColoniesGreen;
xyCorrectionG =  coefsG(1) ...
    + coefsG(2).*xG' + coefsG(3).*yG' ...
    + coefsG(4).*canaryGreen.*canaryMaskRed ...
    + coefsG(5).*distToCenter ...
    + coefsG(6).*xG'.*yG' ...
    + coefsG(7).*xG'.*canaryGreen.*canaryMaskRed ...
    + coefsG(8).*yG'.*canaryGreen.*canaryMaskRed ...
    + coefsG(9).*distToCenter.*canaryGreen.*canaryMaskRed;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize the correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    figure; subplot(2,2,1); imshow(uint16(canaryRedRaw)); 
            xlabel('x'), ylabel('y'), 
            title('Original Canary Red')
            
            subplot(2,2,2); imshow(uint16(xyCorrectionR)); 
            xlabel('x'), ylabel('y'), 
            title('Red With Correction')
    
            subplot(2,2,3); imshow(uint16(canaryGreenRaw)); 
            xlabel('x'), ylabel('y'), 
            title('Original Canary Green') 
            
            subplot(2,2,2); imshow(uint16(xyCorrectionG)); 
            xlabel('x'), ylabel('y'), 
            title('Red With Correction')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the corrected image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFlag
    hRed = uint16(xyCorrectionR);
    pathRed   = strcat(expPath, strcat('/red-corrected/', strcat(char(filenameRed), '.mat')));
    save(pathRed, 'hRed')
    
    hGreen = uint16(xyCorrectionG);
    pathGreen = strcat(expPath, strcat('/green-corrected/', strcat(char(filenameGreen), '.mat')));
    save(pathGreen, 'hGreen')
end

end