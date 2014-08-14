clear all

% For other waveforms, just change the name of xls file, specify the
% axis([]) values for the plots - they would be common for figure 1 (reduced samples and predicted values), and
% then common for figure 2 (error values)
ecgOutpPredbyNN = xlsread('mgecgoutpRED100.xls'); % Predicted output for ART for 598 samples (~600)

ts1 = ecgOutpPredbyNN;
ts1(2:2:end) = [];       % remove alternate samples
samples = 1:length(ts1); % Get the #samples
samples = samples * 2;               % * 2 just to help the plot become similar to predicted one with double the samples by assuming
                                     % that the sample values correspond to
                                     % alternate samples. We need to do
                                     % this to create spaces for inserting
                                     % the predicted values as if they were the missing
                                     % samples
subplot(4,2,1)
plot(samples, ts1, 'k')
title('Original Series with half the samples removed')
axis([0,250,-0.6,1])

yi_ts1 = interp1(samples./2, ts1, 1.5:1:length(ts1), 'nearest'); % Samples need to be set back to actual #samples; 1.5:1: ensures we predict the right middle values by interpolation
yi_ts1(length(yi_ts1)+1) = yi_ts1(length(yi_ts1));
ts1x2 = ts1;                    %make a copy of time series 1 just to show it as a different variable with twice as many points
for j = length(yi_ts1):-1:1    
    ts1x2(2*j)=ts1x2(j);    % each sample copied to twice and twice minus one to create spaces by means of duplicates
    ts1x2(2*j-1)=ts1x2(j);
end
%ts1;
for k=1:length(yi_ts1)
    ts1x2(2*k) = yi_ts1(k);               % insert predicted values in the duplicate valued spaces just created
end
%ts2x2(length(ts2x2)) = [];  % remove the last sample with duplicate values after inserting and completng the time series
samples1 = 1:length(ts1x2);           % Get the new value of samples ~double
subplot(4,2,2)
plot(samples1, ts1x2, 'r')
title('Result of prediction on original Series with half the samples removed')
axis([0,250,-0.6,1])

%disp('#original samples in ts2')
length(ecgOutpPredbyNN);
%disp('#samples in ts2 after sample reduction')
length(ts1);
%disp('#samples in yi')
length(yi_ts1);
%disp('#samples in ts1x2')
length(ts1x2);

err_ts1 = gsubtract(ts1x2, ecgOutpPredbyNN);
ersq_ts1 = err_ts1.*err_ts1;
ersqsum_ts1 = sum(ersq_ts1);
mserr_ts1 = sqrt(ersqsum_ts1/size(err_ts1(:),1));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ts2 = ecgOutpPredbyNN;

originalSamples = ts2;
filledOS = originalSamples; %preserve a copy of the original series

reducedSamples = originalSamples(1:3:end) ;
copyOfOrigReducedSamples = reducedSamples; %copy needed as reducedSamples will change

% If the last element of the series is NOT in the reduced set, we need to
% fill it up with the last sample value
if length(originalSamples) > (1+(length(reducedSamples)-1)*3)
    diff = length(originalSamples) - (1+(length(reducedSamples)-1)*3);
    if diff == 1   % We need to fill up two additional samples at the end with the copy of the last sample in the original series
        for k = 1:diff
            originalSamples(length(originalSamples)+1) = originalSamples(length(originalSamples));
        end
    end
    if diff == 2    % We need to fill up just one additional sample
        originalSamples(length(originalSamples)+1) = originalSamples(length(originalSamples));
    end
end
% recreate the reduced Sample set to include the additional new last
% element in the series (just to avoid a NaN as the last element resulting in subsequent interpolation)  
reducedSamples = originalSamples(1:3:end);
% assign these to every third sample
totalSamples = 1:length(originalSamples);
sampleInstants = totalSamples(1:3:end); % We want to map the samples with x = 1,4,7...

subplot(4,2,3)
plot(sampleInstants, reducedSamples, 'k')
title('Original Series with two/third of the samples removed')
axis([0,250,-0.6,1])

%interpolating for positions 2,5,8.....
yi_2 = interp1(sampleInstants, reducedSamples, 2:3:length(originalSamples), 'nearest');

% The last element in this series containing the result of interpolation has to be removed if we were required to extend the
% original series earlier
if length(originalSamples) > (1+(length(copyOfOrigReducedSamples)-1)*3)
    yi_2(length(yi_2)) = [];
end

%interpolating for positions 2,5,8.....
yi_3 = interp1(sampleInstants, reducedSamples, 3:3:length(originalSamples), 'nearest');
% The last element in this series containing the result of interpolation has to be removed if we were required to extend the
% original series earlier

if length(originalSamples) > (1+(length(copyOfOrigReducedSamples)-1)*3)
    yi_3(length(yi_3)) = [];
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
filledOS;

% Now populating the series with predictions on missing 2nd and 3rd samples 
for count = 1:length(yi_2)
    filledOS((1+(count-1)*3)+1) = yi_2(count);
end

for count = 1:length(yi_3)
    filledOS((1+(count-1)*3)+2) = yi_3(count);
end

filledOS;
subplot(4,2,4)
plot(1:length(filledOS),filledOS, 'm')
title('Result of prediction on original Series with two/third of the samples removed')
axis([0,250,-0.6,1])

err_ts2 = gsubtract(filledOS, ecgOutpPredbyNN);
ersq_ts2 = err_ts2.*err_ts2;
ersqsum_ts2 = sum(ersq_ts2);
mserr_ts2 = sqrt(ersqsum_ts2/size(err_ts2(:),1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

originalSamples = ecgOutpPredbyNN;
filledOS1 = originalSamples; %preserve a copy of the original series

reducedSamples = originalSamples(1:4:end); 
copyOfOrigReducedSamples = reducedSamples; %copy needed as reducedSamples will change


% If the last element of the series is NOT in the reduced set, we need to
% fill it up with the last sample value
if length(originalSamples) > (1+(length(reducedSamples)-1)*4)
    diff = length(originalSamples) - (1+(length(reducedSamples)-1)*4);
    if diff == 1 | diff == 2   % We need to fill up two additional samples at the end with the copy of the last sample in the original series
        for k = 1:diff
            originalSamples(length(originalSamples)+1) = originalSamples(length(originalSamples));
        end
    end
    if diff == 3    % We need to fill up just one additional sample
        originalSamples(length(originalSamples)+1) = originalSamples(length(originalSamples));
    end
end
% recreate the reduced Sample set to include the additional new last
% element in the series (just to avoid a NaN as the last element resulting in subsequent interpolation)  
reducedSamples = originalSamples(1:4:end);
% assign these to every third sample
totalSamples = 1:length(originalSamples);
sampleInstants = totalSamples(1:4:end); % We want to map the samples with x = 1,4,7...

subplot(4,2,5)
plot(sampleInstants, reducedSamples, 'k')
title('Original Series with three/fourth of the samples removed')
axis([0,250,-0.6,1])

%interpolating for positions 2,6,10.....
yi_2 = interp1(sampleInstants, reducedSamples, 2:4:length(originalSamples),'nearest');

% The last element in this series containing the result of interpolation has to be removed if we were required to extend the
% original series earlier
if length(originalSamples) > (1+(length(copyOfOrigReducedSamples)-1)*4)
    yi_2(length(yi_2)) = [];
end

%interpolating for positions 3,7,11.....
yi_3 = interp1(sampleInstants, reducedSamples, 3:4:length(originalSamples), 'nearest');
% The last element in this series containing the result of interpolation has to be removed if we were required to extend the
% original series earlier

if length(originalSamples) > (1+(length(copyOfOrigReducedSamples)-1)*4)
    yi_3(length(yi_3)) = [];
end

%interpolating for positions 4,8,12.....
yi_4 = interp1(sampleInstants, reducedSamples, 4:4:length(originalSamples), 'nearest');
% The last element in this series containing the result of interpolation has to be removed if we were required to extend the
% original series earlier

if length(originalSamples) > (1+(length(copyOfOrigReducedSamples)-1)*4)
    yi_4(length(yi_4)) = [];
end


% The following three for loops replace the alternate 2nd and 3rd samples
% from the copy of the original series with zeroes before populating them
% with the predicted values
for count = 1:length(yi_2)
    filledOS1((1+(count-1)*4)+1) = 0;
end

for count = 1:length(yi_3)
    filledOS1((1+(count-1)*4)+2) = 0;
end

for count = 1:length(yi_4)
    filledOS1((1+(count-1)*4)+3) = 0;
end

% check if the above two for loops worked as expected
filledOS1;

% Now populating the series with predictions on missing 2nd and 3rd samples 
for count = 1:length(yi_2)
    filledOS1((1+(count-1)*4)+1) = yi_2(count);
end

for count = 1:length(yi_3)
    filledOS1((1+(count-1)*4)+2) = yi_3(count);
end

for count = 1:length(yi_4)
    filledOS1((1+(count-1)*4)+3) = yi_4(count);
end

filledOS1;

subplot(4,2,6)
plot(1:length(filledOS),filledOS, 'c')
title('Result of prediction on original Series with three/fourth of the samples removed')
axis([0,250,-0.6,1])

err_ts3 = gsubtract(filledOS1, ecgOutpPredbyNN);
ersq_ts3 = err_ts3.*err_ts3;
ersqsum_ts3 = sum(ersq_ts3);
mserr_ts3 = sqrt(ersqsum_ts3/size(err_ts3(:),1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

originalSamples = ecgOutpPredbyNN;
filledOS2 = originalSamples; %preserve a copy of the original series

reducedSamples = originalSamples(1:5:end); 
copyOfOrigReducedSamples = reducedSamples; %copy needed as reducedSamples will change


% If the last element of the series is NOT in the reduced set, we need to
% fill it up with the last sample value
if length(originalSamples) > (1+(length(reducedSamples)-1)*5)
    diff = length(originalSamples) - (1+(length(reducedSamples)-1)*5);
    if diff == 1 | diff == 2 | diff == 3   % We need to fill up two additional samples at the end with the copy of the last sample in the original series
        for k = 1:diff
            originalSamples(length(originalSamples)+1) = originalSamples(length(originalSamples));
        end
    end
    if diff == 4    % We need to fill up just one additional sample
        originalSamples(length(originalSamples)+1) = originalSamples(length(originalSamples));
    end
end
% recreate the reduced Sample set to include the additional new last
% element in the series (just to avoid a NaN as the last element resulting in subsequent interpolation)  
reducedSamples = originalSamples(1:5:end);
% assign these to every third sample
totalSamples = 1:length(originalSamples);
sampleInstants = totalSamples(1:5:end); % We want to map the samples with x = 1,4,7...

subplot(4,2,7)
plot(sampleInstants, reducedSamples, 'k')
title('Original Series with four/fifth of the samples removed')
axis([0,250,-0.6,1])

%interpolating for positions 2,7,12.....
yi_2 = interp1(sampleInstants, reducedSamples, 2:5:length(originalSamples), 'nearest');

% The last element in this series containing the result of interpolation has to be removed if we were required to extend the
% original series earlier
if length(originalSamples) > (1+(length(copyOfOrigReducedSamples)-1)*5)
    yi_2(length(yi_2)) = ceil(yi_2(length(yi_2)));
end

%interpolating for positions 3,8,13.....
yi_3 = interp1(sampleInstants, reducedSamples, 3:5:length(originalSamples), 'nearest');
% The last element in this series containing the result of interpolation has to be removed if we were required to extend the
% original series earlier

if length(originalSamples) > (1+(length(copyOfOrigReducedSamples)-1)*5)
    yi_3(length(yi_3)) = ceil(yi_3(length(yi_3)));
end

%interpolating for positions 4,9,14.....
yi_4 = interp1(sampleInstants, reducedSamples, 4:5:length(originalSamples), 'nearest');
% The last element in this series containing the result of interpolation has to be removed if we were required to extend the
% original series earlier

if length(originalSamples) > (1+(length(copyOfOrigReducedSamples)-1)*5)
    yi_4(length(yi_4)) = ceil(yi_4(length(yi_4)));
end

%interpolating for positions 5,10,15.....
yi_5 = interp1(sampleInstants, reducedSamples, 5:5:length(originalSamples), 'nearest');
% The last element in this series containing the result of interpolation has to be removed if we were required to extend the
% original series earlier

if length(originalSamples) > (1+(length(copyOfOrigReducedSamples)-1)*5)
    yi_5(length(yi_5)) = ceil(yi_5(length(yi_5)));
end


% The following four for loops replace the alternate 2nd and 3rd samples
% from the copy of the original series with zeroes before populating them
% with the predicted values
for count = 1:length(yi_2)
    filledOS2((1+(count-1)*5)+1) = 0;
end

for count = 1:length(yi_3)
    filledOS2((1+(count-1)*5)+2) = 0;
end

for count = 1:length(yi_4)
    filledOS2((1+(count-1)*5)+3) = 0;
end

for count = 1:length(yi_5)
    filledOS2((1+(count-1)*5)+4) = 0;
end

% check if the above two for loops worked as expected
filledOS2;

% Now populating the series with predictions on missing 2nd and 3rd samples 
for count = 1:length(yi_2)
    filledOS2((1+(count-1)*5)+1) = yi_2(count);
end

for count = 1:length(yi_3)
    filledOS2((1+(count-1)*5)+2) = yi_3(count);
end

for count = 1:length(yi_4)
    filledOS2((1+(count-1)*5)+3) = yi_4(count);
end

for count = 1:length(yi_5)
    filledOS2((1+(count-1)*5)+4) = yi_5(count);
end


filledOS2;

subplot(4,2,8)
plot(1:length(filledOS2),filledOS2, 'b')
title('Result of prediction on original Series with four/fifth of the samples removed')
axis([0,250,-0.6,1])

err_ts4 = gsubtract(filledOS2(1:length(ecgOutpPredbyNN)), ecgOutpPredbyNN);

% Replace NaN values in the series with maximum possible error
for i = 1:length(err_ts4)
    if isnan(err_ts4(i))
        err_ts4(i) = max(err_ts4);
    end
end
 

ersq_ts4 = err_ts4.*err_ts4;
ersqsum_ts4 = sum(ersq_ts4) + 2*(max(err_ts4)*max(err_ts4));
mserr_ts4 = sqrt(ersqsum_ts4/size(err_ts4(:),1));

s_No = 1:4;
no_of_samples = ['    Halved' '     1/3rd' '     1/4th' '     1/5th']
mse = [mserr_ts1 mserr_ts2 mserr_ts3 mserr_ts4]

figure(2)
subplot(2,2,1)
plot(1:length(err_ts1),err_ts1, 'r')
title('Error plot for prediction on original Series with half of the samples removed')
axis([0,250,-0.3,0.3])

subplot(2,2,2)
plot(1:length(err_ts2),err_ts2, 'm')
title('Error plot for prediction on original Series with two/third of the samples removed')
axis([0,250,-0.3,0.3])

subplot(2,2,3)
plot(1:length(err_ts3),err_ts3, 'c')
title('Error plot for prediction on original Series with three/fourth of the samples removed')
axis([0,250,-0.3,0.3])

subplot(2,2,4)
plot(1:length(err_ts4),err_ts4, 'b')
title('Error plot for prediction on original Series with four/fifth of the samples removed')
axis([0,250,-0.3,0.3])

disp('--------------------')
disp('Signal = ECG')
disp('--------------------')
disp('Signal minimum = ') 
disp(min(ecgOutpPredbyNN))
disp('Signal maximum = ') 
disp(max(ecgOutpPredbyNN))
disp('Error minimum = for the four cases') 
disp(min(err_ts1)) 
disp(min(err_ts2)) 
disp(min(err_ts3)) 
disp(min(err_ts4))
disp('Error maximum = for the four cases') 
disp(max(err_ts1)) 
disp(max(err_ts2)) 
disp(max(err_ts3)) 
disp(max(err_ts4))
disp('--------------------')

