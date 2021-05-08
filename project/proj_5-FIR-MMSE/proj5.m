clear; close all; clc;

% c[n]
c = [1 .2 .4];

% simulation for N = 4, 6, 10
N = [4,6,10];

% array to store MSE of each N
MSE = zeros(1,length(N));

% signal length
len_N = 1000;

% random signal (+/-1)
s = randi(2, 1, len_N);
s(s == 2) = -1;

[Rss, lags] = xcorr(s); % index 1000 is time delay = 0;

% output of 1st filter c[n]
y = filter(c, 1, s);

% standard deviation of noise
sd = 0.5;

% y + d[n] -> input of 2nd filter h[n]
r = y + normrnd(0, sd, 1, length(y));

% simulate for each N
for m = 1:length(N)
    % Rrr -> autocorrelation
    Rrr = xcorr(r); 
    Rrr_mid = (length(Rrr)+1)/2; 
    % Rsr
    Rsr = xcorr(s,r); 
    Rsr_mid = (length(Rsr)+1)/2; 
    Left = zeros(N(m));
    for i = 1:N(m)
        Left(i,:) = transpose(Rrr(Rrr_mid - i+1 : Rrr_mid - i+N(m))); 
    end
    right = transpose(Rsr(Rsr_mid : Rsr_mid + N(m)-1)); 
    
    % slove for h[n]
    h = Left\right;
    
    % output of 2nd filter h[n]
    s_hat = filter(transpose(h), 1, r);
    
    % MSE
    MSE(m) = sum((s_hat-s).^2)/len_N; 
end
table(MSE(1), MSE(2), MSE(3), 'VariableNames', ["N = 4", "N = 6", "N = 10"],...
 'RowNames', "MSE")