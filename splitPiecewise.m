function [lowerConstraints, upperConstraints] = splitPiecewise(oldConstraints, splitPoint)

lowerConstraints = {};
upperConstraints = {};


for constraint = 1:2
    
    for month=1:size(oldConstraints{constraint}(:,:,:),3)
        
        belowSplit = oldConstraints{constraint}(1,:,month) <= splitPoint(1);
        aboveSplit = ~belowSplit;
        
        lower = oldConstraints{constraint}(:,belowSplit, month);
        upper = oldConstraints{constraint}(:,aboveSplit, month);
        
        lowerConstraints{constraint}(:,:,month) = [lower splitPoint'];
        upperConstraints{constraint}(:,:,month) = [splitPoint' upper];
    end
end


end