% Finds plate edge and removes the outer radius of the plate up to
% the designated disk size
%
% Required inputs: Binary plate image (?),
%                  disksize [For erosion radius]
% Optional inputs: Visualization flag
function erosion = erodePlate(preEroded, plateCenter, plateRadius, maxRadiusPercent, visualFlag)

if nargin < 4
    maxRadiusPercent = 0.85;
else
    if nargin < 5
        visualFlag = 0;
    end
end

mask = zeros(size(preEroded));

minDistance = plateRadius * maxRadiusPercent;
z = mask;

for i = 1:size(preEroded,1)
    for j = 1:size(preEroded,2)
        z(i,j) = sqrt((plateCenter(1)-j)^2 + (plateCenter(2)-i)^2);
        if z(i,j) > minDistance
            mask(i,j) = 0;
        else
            mask(i,j) = 1;
        end
    end
end

erosion = logical(mask);%preEroded.*logical(mask);


if visualFlag
    figure
    title('Mask!')
    imshow(erosion)
end

end