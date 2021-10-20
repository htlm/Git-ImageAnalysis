%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Function to run pkfnd with on a specific frame, threshold and sz
%%%%%% parameter.
%%%%%%
%%%%%% filenameArray - is an array with full file locations to red and green
%%%%%% .mat files for the experiment. 
%%%%%%
%%%%%% paramInput    - is a triplet of (frame #, threshold, sz) parameters.
%%%%%%               - it gets separated into the paramArray structure
%%%%%%
%%%%%% plateParams   - is a struct with two subfields that depend on the
%%%%%%                 dataset
%%%%%%                 plateParams.radius - plate radius
%%%%%%                 plateParams.center - [x, y] coordinates
%%%%%%
%%%%%% Last updated: 08/09/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pks = fnRunPkFnd(filenameArray, expPath, paramInput, plateParams)
imageFlag = 0;
plateErosionParam = 0.9; % Arbitrary
paramArray.frame = paramInput(1);
paramArray.th    = paramInput(2);
paramArray.sz    = paramInput(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load images - only need the red data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filenameR    = strcat(filenameArray{1}, num2str(paramArray.frame), '.mat');
filenameREnd = strcat(filenameArray{1}, num2str(plateParams.numImages), '.mat');

load(filenameREnd);
hRedEnd = hRed;
load(filenameR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a binary mask with an eroded the plate edge to minimize artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
erodedMask = erodePlate(imbinarize(hRed), plateParams.center, plateParams.radius, ...
                        plateErosionParam, 0);

if imageFlag
    figure(1)
    subplot(1,2,1), imshow(imadjust(hRed))    
    subplot(1,2,2), imshow(imadjust(hRed.*uint16(erodedMask)))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pks = pkfnd(hRed.*uint16(erodedMask), paramArray.th, paramArray.sz);

h = figure(2);
imshow(hRedEnd); hold on
if ~isempty(pks)
    scatter(pks(:,1), pks(:,2))
    title('')
    csvwrite(strcat(expPath, '/pks_', num2str(paramArray.frame), '-', ...
             num2str(paramArray.th), '-', num2str(paramArray.sz), ...
             '.csv'),...
             pks);
else
    title('no peaks')
end
saveas(h, strcat(expPath, '/pks_', num2str(paramArray.frame), '-', ...
                 num2str(paramArray.th), '-', num2str(paramArray.sz), ...
                 '.png'));


end