% After running trainctl01.m, run the following code for the plot
plot(tamp,errors1,'g')
% Run initially with the following statement commented
% see the error, then decide on the latter two values
axis([0,3600,-0.2,0.3])
hold on;
errPID = simplePID(:)-snI(:);
plot(samp,errPID,'r')
xlabel('Samples')
ylabel('Error')
%Change the name of the plot below
title('PAP: Error plot of PID v/s Non-linear Regression based prediction')
legend('NLR','PID')
