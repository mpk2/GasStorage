function [lowerConstraints, upperConstraints] = splitPiecewise(oldConstraints, splitPoint, invalidMonthIndex)

lowerConstraints = {};
upperConstraints = {};


for constraint = 1:2
    
    for month=1:length({oldConstraints{constraint,:}})
        
        if(month == invalidMonthIndex) 
            belowSplit = oldConstraints{constraint,month}(1,:) <= splitPoint(1);
            aboveSplit = ~belowSplit;

            lower = oldConstraints{constraint, month}(:,belowSplit);
            upper = oldConstraints{constraint, month}(:,aboveSplit);

            lowerConstraints{constraint, month} = [lower splitPoint'];
            upperConstraints{constraint, month} = [splitPoint' upper];
            
        else
            lowerConstraints{constraint, month} = oldConstraints{constraint, month};
            upperConstraints{constraint, month} = oldConstraints{constraint, month};
        end
    end
end


end