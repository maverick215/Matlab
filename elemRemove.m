%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program works
% It removes every 2nd and 3rd sample from the original series and copies
% the rest of the samples into a reduced samples vector. The missing values
% are then predicted and the original series is reconstructed by filling up
% with the predicted values
% Remove every 2nd and 3rd value 
originalSamples = 1:12
filledOS = originalSamples; %preserve a copy of the original series

reducedSamples = originalSamples(1:3:end) 
copyOfOrigReducedSamples = reducedSamples; %copy needed as reducedSamples will change

% If the last element of the series is NOT in the reduced set, we need to
% fill it up with the last sample value
if length(originalSamples) > (1+(length(reducedSamples)-1)*3)
    diff = length(originalSamples) - (1+(length(reducedSamples)-1)*3)
    if diff == 1   % We need to fill up two additional samples at the end with the copy of the last sample in teh original series
        for k = 1:diff
            originalSamples(length(originalSamples)+1) = originalSamples(length(originalSamples))
        end
    end
    if diff == 2    % We need to fill up just one additional sample
        originalSamples(length(originalSamples)+1) = originalSamples(length(originalSamples))
    end
end
% recreate the reduced Sample set to include the additional new last
% element in the series (just to avoid a NaN as the last element resulting in subsequent interpolation)  
reducedSamples = originalSamples(1:3:end)
% assign these to every third sample
totalSamples = 1:length(originalSamples)
sampleInstants = totalSamples(1:3:end) % We want to map the samples with x = 1,4,7...

%interpolating for positions 2,5,8.....
yi_2 = interp1(sampleInstants, reducedSamples, 2:3:length(originalSamples))

% The last element in this series containing the result of interpolation has to be removed if we were required to extend the
% original series earlier
if length(originalSamples) > (1+(length(copyOfOrigReducedSamples)-1)*3)
    yi_2(length(yi_2)) = []
end

%interpolating for positions 2,5,8.....
yi_3 = interp1(sampleInstants, reducedSamples, 3:3:length(originalSamples))
% The last element in this series containing the result of interpolation has to be removed if we were required to extend the
% original series earlier

if length(originalSamples) > (1+(length(copyOfOrigReducedSamples)-1)*3)
    yi_3(length(yi_3)) = []
end


% The following two for loops replace the alternate 2nd and 3rd samples
% from the copy of the original series with zeroes before populating them
% with the predicted values
for count = 1:length(yi_2)
    filledOS((1+(count-1)*3)+1) = 0;
end

for count = 1:length(yi_3)
    filledOS((1+(count-1)*3)+2) = 0;
end

% check if the above two for loops worked as expected
filledOS

% Now populating the series with predictions on missing 2nd and 3rd samples 
for count = 1:length(yi_2)
    filledOS((1+(count-1)*3)+1) = yi_2(count);
end

for count = 1:length(yi_3)
    filledOS((1+(count-1)*3)+2) = yi_3(count);
end

filledOS
