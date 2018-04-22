%data is a data matrix with each rows for each subject and columns for that
%subjects score on each variable.

alpha = .05;
n = size(data,1);
p = size(data,2);

for i_var=1:(p-1)
    tran(:,i_var) = data(:,i_var)-data(:,p);
end

X_bar_D = mean(tran,1);
mu_D = zeros(1,(p-1));
S_D = cov(tran);


T_square = n* ((X_bar_D - mu_D) * inv(S_D) * (X_bar_D - mu_D)')

T_square_crit = ((n-1)*(p-1) / (n-p+1))*(finv(1-alpha,p-1,n-p-1))







