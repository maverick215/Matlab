% Reads the data from TestDataCalVal.xlsx file [128x13] and plots it in line connections

testmat = xlsread('TestDataCalVal');
x = [1:128];
orig = testmat(:,1);
plot(x,orig);
hold on;
u = testmat(:,2);
z = testmat(:,5);
w = testmat(:,8);
v = testmat(:,11);
plot(x,u, 'y')
plot(x,z, 'r');
plot(x,w, 'k');
plot(x,v, 'g');
ylabel('EEG')
xlabel('Sample #')
legend('Orig', 'Linear', 'PID', 'MVAVG', 'Past')
axis([2,128,1750,2150])
title('Comparison of Original Series with predicted values from four algorithms') 
