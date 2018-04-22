function [F r_diff] = circle_test(cos_out,sin_out,cos_bar,sin_bar)
%det-undet (2-1)
N = size(cos_out,1);

X_diff = cos_out(:,2)-cos_out(:,1);
Y_diff = sin_out(:,2)-sin_out(:,1);
r_diff = sqrt ( (X_diff.^2) + (Y_diff.^2) );



X_bar_diff = cos_bar(1,2)-cos_bar(1,1);
Y_bar_diff = sin_bar(1,2)-sin_bar(1,1);
SS_X = (sum(X_diff.^2)) - ((sum(X_diff).^2)/N);
SS_Y = (sum(Y_diff.^2)) - ((sum(Y_diff).^2)/N);
SS_XY = sum(X_diff.*Y_diff) - ((sum(X_diff).*sum(Y_diff))/N);

num = ( (X_bar_diff^2)*SS_Y) - (2*( X_bar_diff*Y_bar_diff*SS_XY)) + ((Y_bar_diff^2)*SS_X);
denom = (SS_X*SS_Y)-(SS_XY^2);
constant = (N*(N-2))/2;


F = constant * (num/denom);
% 
% 
% 








% cos_A = X_diff./r_diff;
% sin_A = Y_diff./r_diff;
% 
% ranker(:,1) = r_diff;
% ranker(:,2) = cos_A;
% ranker(:,3) = sin_A;
% ranker = sortrows(ranker,1);
% for x=1:size(ranker,1)
%     ranker(x,4) = x;
% end
% 
% super_X = sum(ranker(:,4).*cos(ranker(:,3))) / N;
% super_Y = sum(ranker(:,4).*sin(ranker(:,2))) / N;
% super_r = sqrt (( (super_X^2) + (super_Y^2)) / N)
% 
% A_diff = atand(Y_diff./X_diff);
% 
% grand_diff = circle_grand_mean(A_diff,r_diff)

%-----------------------





% X_bar_diff = X_bar(1,2)-X_bar(1,1);
% Y_bar_diff = Y_bar(1,1)-X_bar(1,1);
% 
% if Y_bar_diff > 0 & X_bar_diff > 0
%     circ_mean_diff = atand(Y_bar_diff/X_bar_diff);
% elseif X_bar_diff < 0
%     circ_mean_diff = atand(Y_bar_diff/X_bar_diff) + 180;
% elseif Y_bar_diff < 0 & X_bar_diff > 0 
%     circ_mean_diff = atand(Y_bar_diff/X_bar_diff) + 360;
% end
% 
% K = size(X_out,1);
% 
% SS_X = sum((X_diff.^2)) - ((sum(X_diff).^2)/K);
% SS_Y = sum((Y_diff.^2)) - ((sum(Y_diff).^2)/K);
% SS_XY = sum(X_diff.*Y_diff) - (sum(X_diff).*sum(Y_diff));
% 
% num = ( (X_bar_diff^2)*SS_Y) - (2*( X_bar_diff*Y_bar_diff*SS_XY)) + ((Y_bar_diff^2)*SS_X);
% denom = (SS_X*SS_Y)-(SS_XY^2);
% constant = (K*(K-2))/2;
% 
% 
% F = constant * (num/denom);



