%Version 3 with LinCheck Adaptive Sampling Algorithm

% Reads the data from TestDataCalVal.xlsx file [128x13] and plots it in line connections

clear all
simplenarxInputs = xlsread('apneaECG400.xls');
simplenarxInputs = medfilt1(simplenarxInputs)';
%simplenarxInputs = simplenarxInputs(1:3:size(simplenarxInputs, 2));

%make a copy of matrix
snI = simplenarxInputs;

%LinCheck Algorithm
for i = 1 : (length(snI) - 2)
    if (snI(i+2) - snI(i+1)) == (snI(i+1) - snI(i))
        simplenarxInputs(i+2) = 0;
    end
end

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
cell2csv('apn400Alterr.csv', errors);
cell2csv('apn400Altoutp.csv', y);
perf = mse(net, Ts, y)
