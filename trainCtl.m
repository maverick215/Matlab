% Reads the data from .xlsx file and plots it in line connections

clear all
% Change the names of the below sheets
simplenarxInputs = xlsread('MGHMFsamples.xls','ecgO');
simplePID = xlsread('MGHMFsamples.xls','ecgP');

%if there is a preparets error, we need to put the ' at the end to
%transpose
simplenarxInputs = medfilt1(simplenarxInputs)';
% %simplenarxInputs = simplenarxInputs(1:3:size(simplenarxInputs, 2));

simplenarxTargets = simplenarxInputs;
simplenarxInputs = num2cell(simplenarxInputs);
simplenarxTargets = simplenarxInputs;
T = simplenarxTargets;
X = simplenarxTargets;
net = narxnet(1:2, 1:2, 10);
[Xs, Xi, Ai, Ts] = preparets(net, X, {}, T);
net = trainlm(net,Xs,Ts,Xi,Ai);
y = net(Xs,Xi,Ai);
errors = gsubtract(Ts,y);
% change the names of below files
cell2csv('mgecgErr.csv', errors);
cell2csv('mgecgoutp.csv', y);
perf = mse(net, Ts, y)
errors1 = cell2mat(errors);
ersq = errors1.*errors1;
ersqsum = sum(ersq)
ermsq = sqrt(ersqsum/size(errors1(:),1))
simplenarxInputs = cell2mat(simplenarxInputs);
samp = 1:length(simplenarxInputs);
plot(samp,simplenarxInputs,'b');
hold on;
plot(samp, simplePID,'r');
y = cell2mat(y);
tamp = 1:length(y);
plot(tamp, y, 'g')
% Change the latter two values
axis([0,3600,-0.7,1])
xlabel('Samples')
ylabel('Predicted Values')
% Change the name of the parameter below 
title('ECG: PID v/s Non-linear Regression based prediction')
legend('Inputs','PID','NLR')

% Code to plot error
% After running trainctl01.m, run the following code for the plot
% plot(tamp,errors1,'g')
% axis([0,3600,-2,2])
% hold on;
% errPID = simplePID(:)-simplenarxInputs(:);
% plot(samp,errPID,'r')
% xlabel('Samples')
% ylabel('Error')
% title('CVP: Error plot of PID v/s Non-linear Regression based prediction')
% legend('PID','NLR')

% Next we go to find how many samples would not be transmitted
% and would be saved because of prediction

%make a copy of matrix
snI = simplenarxInputs;

%LinCheck Algorithm
for i = 1 : (length(snI) - 2)
    if (snI(i+2) - snI(i+1)) <= 0.1 %change this value
        simplenarxInputs(i+2) = 0;
    end
end
sumzeros = sum(simplenarxInputs(:) == 0)

% Be sure to save the time series response of NN algorithm, prediction
% comparison graph, error comparison graph, 
% and values of --- #iterations, perf, ersqsum,
% ermsq, sumzeros
