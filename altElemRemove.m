% % How to remove alternate elements from an array
% m = [1 2 3 4 5 6 7 8];
% m(2:2:end) = []


% This programs tries to find the alternate missing samples
% and see how close the resulting prediction is to the original time series
% 
% The following commented piece is the test code
% to elongate a vector i and stuff another vector j in between
% i = [1 3 5 7 9];
% j = [2 4 6 8];
% 
% for k=length(i):-1:1
%     i(2*k)=i(k);
%     i(2*k-1)=i(k);
% end
% 
% for l=1:length(j)
%     i(2*l)=j(l);
% end
% 
% i(length(i)) = [];

