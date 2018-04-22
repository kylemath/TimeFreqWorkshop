function [circ_mean,range_bar,X_bar,Y_bar,cos_bar,sin_bar] = circle_grand_mean(data,range_data)

%finds the grand mean of the set of circular means in one row of data,
%considering their range data in the range_data

N = size(data,1);

Y_bar = sum(range_data(:,1).* sind(data(:,1)) )/ N; %
X_bar = sum(range_data(:,1).* cosd(data(:,1)) )/ N; %


range_bar = sqrt((X_bar^2)+(Y_bar^2));

if Y_bar > 0 & X_bar > 0
    circ_mean = atand(Y_bar/X_bar);
elseif X_bar < 0
    circ_mean = atand(Y_bar/X_bar) + 180;
elseif Y_bar < 0 & X_bar > 0 
    circ_mean = atand(Y_bar/X_bar) + 360;
end

cos_bar = X_bar/range_bar;
sin_bar = Y_bar/range_bar;