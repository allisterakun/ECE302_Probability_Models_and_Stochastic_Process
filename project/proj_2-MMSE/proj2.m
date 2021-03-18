% ECE 302 - Probability Models and Stochastic Process
% Project 2 - MMSE
% Min (Ella) Cheng, Amy Leong, Allister Liu

%%
clc; clear; close all;
% MMSE Estimation
% Due March 18, 2021 8:00 AM
% 
% Instructions:
% 	MATLAB exercise
% 	Estimation techniques
% 
% Overview: 
%     In this exercise, you will construct several estimators and compare
% the results. You will implement Bayesian and linear MMSE estimators.
% 
% 	Scenario 1:
% 		Implement the Bayes MMSE and Linear MMSE estimators from examples 
%   8.5 and 8.6. Simulate this system by random draws of Y and W, and then 
%   estimating Y from the observations X = Y + W. Verify that your 
%   simulation is correct by comparing theoretical and empirical values of
%   the MSE. Report your results in a table.
% 
% 	Scenario 2: 
% 		Implement the linear estimator for multiple noisy observations, 
%   similar to example 8.8 from the notes. Extend this example so that it 
%   works for an arbitrary number of observations. Use Gaussian random 
%   variables for Y and R. Set μy = 1. Experiment with a few different 
%   variances for both Y and R. On one plot, show the mean squared error 
%   of your simulation compared to the theoretical values for at least 2 
%   different pairs of variances.
%% Scenario 1
clc; clear;

N = 1e6;
Y = 2*rand(1,N)-1;
W = 4*rand(1,N)-2;
X = Y + W;
X_old = X;

X(X>-1 & X<1) = 0;
X(X<-1) = 0.5+X(X<-1)*.5;
X(X>1) = -0.5+X(X>1)*.5;

MSE = mean((Y-X).^2);
bayesMSE = mean(MSE);

linear = X_old/5;
temp = mean((Y-linear).^2);
linearMSE = mean(temp);

T = table([1/4;4/15],[bayesMSE;linearMSE],'VariableNames', ...
          {'Theoretical','Experimental'},'RowNames',{'bayes','linear'});
disp(T)

%% Scenario 2
clc; clear;

theo = zeros(20,4,1);
expr = zeros(20,4,1);
for m = 1:20
    [theo(m,:), expr(m,:)] = scenario2(m);
end
figure
for p = 1:4
   plot(1:1:20, theo(:,p))
   hold on
   scatter(1:1:20, expr(:,p),'o')
   hold on
end
title("MMSE and Number of Observations using \mu_{\it Y} = 1, \mu_{\it R} = 0");
xlabel("Number of Observatons");
ylabel("MMSE");legend(...
    'Theoretical: \sigma_{\it Y}^2 = 0.25, \sigma_{\it R}^2 = 0.25 ', ...
    'Experimental: \sigma_{\it Y}^2 = 0.25,, \sigma_{\it R}^2 = 0.25 ', ...
    'Theoretical: \sigma_{\it Y}^2 = 0.5, \sigma_{\it R}^2 = 0.5', ...
    'Experimental: \sigma_{\it Y}^2 = 0.5, \sigma_{\it R}^2 = 0.5', ...
    'Theoretical: \sigma_{\it Y}^2 = 0.75, \sigma_{\it R}^2 = 0.75', ....
    'Experimental: \sigma_{\it Y}^2 = 0.75, \sigma_{\it R}^2 = 0.75', ...
    'Theoretical: \sigma_{\it Y}^2 = 1, \sigma_{\it R}^2 = 1', ...
    'Experimental: \sigma_{\it Y}^2 = 1, \sigma_{\it R}^2 = 1');

function [theoMMSE, expMMSE] = scenario2(nObs)
    N = 1e6;
    %nObs = 20;
    varY = [0.25,0.5,0.75,1];
    varR = [0.25,0.5,0.75,1];

    % Experimental MMSE 
    expMMSE = zeros(4,1);
    % Theoretical MMSE
    theoMMSE = zeros(4,1);
    for i = 1:4
        theoMMSE(i,:) = (varY(i) * varR(i)) / (nObs * varY(i) + varR(i));
    end

    Y = zeros(4, N, 1);
    R = zeros(4, N, nObs);
    X = zeros(4, N, nObs);
    Yhat = zeros(4, N ,1);

    for i = 1:4
       Y(i,:,:) = normrnd(1, sqrt(varY(i)), [N 1]);
       R(i,:,:) = normrnd(0, sqrt(varR(i)), [N nObs]);
       for j = 1:nObs
          X(i,:,j) = R(i,:,j) + Y(i,:,:); 
       end

       Xtemp = zeros(N, nObs);
       Xtemp(:,:) = X(i,:,:);
       Ytemp = Y(i,:)';

       varYtemp = var(Ytemp);
       muYtemp = mean(Ytemp);
       varianceR = zeros(nObs, 1);
       for k = 1:nObs
          XminusY = Xtemp(:,k)-Ytemp;
          varianceR(k) = var(XminusY); 
       end
       avgVarR = mean(varianceR);

       mms = (1 / (nObs * varYtemp +avgVarR)) * (avgVarR * muYtemp + varYtemp * sum(Xtemp, 2));
       Yhat(i,:,:) = mms;
       expMMSE(i,1) = mean((Ytemp - mms) .^ 2);
    end
end
