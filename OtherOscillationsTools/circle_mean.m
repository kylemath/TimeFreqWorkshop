function [circ_mean,range,X,Y,cos_a,sin_a] = circle_mean(data)

%finds the mean of the set of degrees in one row of data

sumfor_Y = 0;
sumfor_X = 0;

for count = 1:size(data)                                %Don't include the zeros when there arnt 499 trials in this sub/cond
    if data(count,1) == 0
        break
    end
    if data(count,1) ~= 0                       %Sum up all the sin(Y) and cos(X) of the angles
        sumfor_Y = sumfor_Y + sind(data(count,1));
        sumfor_X = sumfor_X + cosd(data(count,1));
    end
end

Y = sumfor_Y / (count-1);           %Divide over the number of trials
X = sumfor_X / (count-1);

range = sqrt((X^2)+(Y^2));          %get the range or length of the mean vector

if Y > 0 & X > 0
    circ_mean = atand(Y/X);
elseif X < 0
    circ_mean = atand(Y/X) + 180;
elseif Y < 0 & X > 0 
    circ_mean = atand(Y/X) + 360;
end

cos_a = X/range;            %get these two other things for later
sin_a = Y/range;

