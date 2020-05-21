% Perform background correction on Canary image dataset
% Typhoon Plate scan from 06282018 used to create correction
%
% Last updated: 08/06/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fnDatasetBkgndCorrectionUpdated(expPath, imageRange)

saveFlag = 1;
pre2017DatasetFlagBlock = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse expPath name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathParse = strsplit(expPath,'/');
expName = char(pathParse(length(pathParse)));
post2017DatasetFlag = 0;

if isempty(expName)
    expName = char(pathParse(length(pathParse) - 1));
end

expIdentifier = char({expName(1:8)});
bkgndCorrYrs = csvread('bkgndCorrYr.csv', 1,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify which background correction to use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty( find(str2double(expIdentifier) == bkgndCorrYrs(:,1)))

    bkgndCorrYr = 2019;
    
else
    I = find(str2double(expIdentifier) == bkgndCorrYrs(:,1));
    if isempty(I)
        I = find(str2double(expIdentifier(2:end)) == bkgndCorrYrs(:,1));
    end
    bkgndCorrYr = bkgndCorrYrs(I,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct filenames 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if str2double(expIdentifier(end-3:end)) > 2017
    post2017DatasetFlag = 1;
end

if bkgndCorrYr > 2017
    post2017DatasetFlag = 1;
end

if post2017DatasetFlag
    filenameArray = {strcat(expName, '_exc_DsRed_em_DsRed_channel2_t');...
    strcat(expName, '_exc_GFP_em_GFP_channel1_t')};
else
    filenameArray = {strcat(expName, '_exc_DsRed_em_DsRed_t');...
    strcat(expName, '_exc_GFP_em_GFP_t')};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Background correct each image with the appropriate correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = imageRange
    
    if i < 100
        imageNum = strcat('00', num2str(i));
        
        if i < 10
            imageNum = strcat('0', imageNum);
        end
    else
        imageNum = strcat('0', num2str(i));
    end
    
    filenameRed   = strcat(filenameArray{1}, imageNum);
    filenameGreen = strcat(filenameArray{2}, imageNum);
    
    fnBackgroundCorrection06202019(filenameRed, filenameGreen, expPath, ...
                                      saveFlag, bkgndCorrYr);
    
%     if pre2017DatasetFlagBlock
%         fnBackgroundCorrectionUpdatedAgain(filenameRed, filenameGreen, expPath, ...
%                                       saveFlag, post2017DatasetFlag);
%     else
%         fnBackgroundCorrection(filenameRed, filenameGreen, expPath, saveFlag);
%     end
    
end

end