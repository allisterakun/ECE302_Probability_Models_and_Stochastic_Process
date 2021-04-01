%Allister Liu, Min (Ella) Cheng, Amy Leong
%Stochastics Project 3

clc; 
clear;
close all;


%% Part 1 
%setup
N = 1000000;
observation = 3:13;
lambda = [1 5 10];
alpha = [1 5 10];

ind=1;
for i = observation
    %calculate the MSE, bias and variance for each distribution and for
    %each # of observation
    
    %MSE is found by taking the mean square of the ML
    %estimator-lambda/alpha
    
    %bias is found by taking the mean of the ML estimator -lambda/alpha
    
    %variance is found by taking the variance of the ML estimator  

    [MSE_exp1(ind),bias_exp1(ind),var_exp1(ind)]= exponential_mse_bias_variance(N,i,lambda(1));
    [MSE_exp2(ind),bias_exp2(ind),var_exp2(ind)]= exponential_mse_bias_variance(N,i,lambda(2));
    [MSE_exp3(ind),bias_exp3(ind),var_exp3(ind)]= exponential_mse_bias_variance(N,i,lambda(3));
    [MSE_rayleigh1(ind),bias_ray1(ind),var_ray1(ind)]=rayleigh_mse_bias_variance(N,i,alpha(1));
    [MSE_rayleigh2(ind),bias_ray2(ind),var_ray2(ind)]=rayleigh_mse_bias_variance(N,i,alpha(2));
    [MSE_rayleigh3(ind),bias_ray3(ind),var_ray3(ind)]=rayleigh_mse_bias_variance(N,i,alpha(3));
    
    ind=ind+1;
end

% plot the MSE for exponential and rayleigh distributions
figure;
subplot(1,2,1);
plot(observation, MSE_exp1, observation, MSE_exp2, observation,MSE_exp3);
title("Exponential MSE");
xlabel("Number of Observations");
ylabel("MSE");
legend("\lambda = " + lambda(1), "\lambda = " + lambda(2),"\lambda = " + lambda(3));

subplot(1,2,2);
plot(observation, MSE_rayleigh1, observation, MSE_rayleigh2, observation, MSE_rayleigh3);
title("Rayleigh MSE");
xlabel("Number of Observations");
ylabel("MSE");
legend("\alpha = " + alpha(1), "\alpha = " + alpha(2),"\alpha = " + alpha(3));

%plot the bias for exponential and rayleigh distributions
figure;
subplot(1,2,1);
plot(observation, bias_exp1, observation, bias_exp2, observation, bias_exp3);
title("Exponential Bias");
xlabel("Number of Observations");
ylabel("Bias");
legend("\lambda = " + lambda(1), "\lambda = " + lambda(2),"\lambda = " + lambda(3));


subplot(1,2,2);
plot(observation, bias_ray1, observation, bias_ray2,observation,bias_ray3);
title("Rayleigh Bias");
xlabel("Number of Observations");
ylabel("Bias");
legend("\alpha = " + alpha(1), "\alpha = " + alpha(2),"\alpha = " + alpha(3));

%plot the variance for exponential and rayleigh distributions
figure;
subplot(1,2,1);
plot(observation, var_exp1, observation, var_exp2, observation, var_exp3);
title("Exponential Variance");
xlabel("Number of Observations");
ylabel("Variance");
legend("\lambda = " + lambda(1), "\lambda = " + lambda(2),"\lambda = " + lambda(3));

subplot(1,2,2);
plot(observation, var_ray1, observation, var_ray2,observation,var_ray3);
title("Rayleigh Variance");
xlabel("Number of Observations");
ylabel("Variance");
legend("\alpha = " + alpha(1), "\alpha = " + alpha(2),"\alpha = " + alpha(3));



%% Part 2
load data.mat;
data1 = data.';     
[~, size] = size(data);    %get the size of data

%calculate the ML estimators
alphaEstimator = sqrt(.5 * mean(data.^2, 2));
lambdaEstimator = size./sum(data,2);

%get the sum of log likelihoods for exponential and rayleigh distribution 
%the pdf of an exponential distribution is given as lambdaEstimator * exp(-lambdaEstimator * data1
%the pdf of a Rayleigh distribution is given as data1/alphaEstimator^2 .*
%exp(-data1.^2/(2*alphaEstimator^2))

exponentialLikelihood=sum(log(lambdaEstimator * exp(-lambdaEstimator * data1)));
rayleighLikelihood= sum(log(data1/alphaEstimator^2 .* exp(-data1.^2/(2*alphaEstimator^2))));

disp("The sum of log likelihoods for exponential distribution is " +exponentialLikelihood)
disp("The sum of log likelihoods for rayleigh distribution is "+rayleighLikelihood)

%The log likelihood value measures how well a model fits. The higher the value, the better
%the fit is. Since the rayleigh distribution has a higher likelihood, the data was most
%likely drawn from a rayleigh distribution. 



function [mse, bias, variance] = rayleigh_mse_bias_variance(N,i,alpha)
    rayleigh = raylrnd(alpha, [N i]);     %getting i samples with random draws
    avg=mean(rayleigh.^2,2);
    alpha2 = sqrt(.5 * avg);          %ML estimator
    mse= mean((alpha - alpha2).^2);   %get MSE
    bias = mean(alpha2) - alpha;   %get bias
    variance = var(alpha2);         %get variance


end


function [mse, bias, variance] = exponential_mse_bias_variance(N,i,lambda)
   
    exponential = exprnd(1/lambda, [N i]); %getting i samples with random draws
    add=sum(exponential,2);
    lambda2 =i ./ add;               %ML estimator
    mse = mean((lambda- lambda2).^2);  %get the MSE 
    bias= mean(lambda2) - lambda;    %get the bias
    variance = var(lambda2);            %get the variance 


end
