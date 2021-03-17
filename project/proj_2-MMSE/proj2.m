
% scenario1
clear;
N = 100000;
Y = -1 + (1-(-1)).*rand(1, N);
W = -2 + (2-(-2)).*rand(1, N);

X = W+Y;
fx = zeros(1, N);
% estimator for conditional probability
bayesEstimator = zeros(1,N);
% sequare error of bayes estimator
bayesSE = zeros(1,N);

for i=1:1:N
   fx(i) = 0;
   if X(i)>=-3 && X(i)<-1
       fx(i) = (3+X(i))/8;
       bayesEstimator(i) = 0.5+0.5*X(i);
   end
   if X(i)>=-1 && X(i)<1
       fx(i) = 0.25;
       bayesEstimator(i) = 0;
   end
   if X(i)>=1 && X(i)<=3
       fx(i) = (3-X(i))/8;
       bayesEstimator(i) = -0.5+0.5*X(i);
   end
   bayesSE = (Y(i)-bayesEstimator(i)).^2;
    
end
% mean square error of Y
bayesMSE = sum(bayesSE)/N;
bMMSE = zeros(1,N);
for i=1:1:N
    bMMSE(i) = bayesMSE * fx(i);
end
% bayesian mean square estimator
MMSE = sum(bMMSE)/N;
%display(MMSE);

% LMSE
meanx = sum(X)/N;
meany = sum(Y)/N;
%covariance based on eqn. 8.19
covXY = sum(X.*Y)/N-meanx*meany;
varX = sum(X.^2)/N-meanx^2;

% linear estimator for y
ylEstimator = zeros(1,N);
% linear estimator square error
ylSE = zeros(1, N);
for i=1:1:N
    %linear estimator based on eqn.8.51
    ylEstimator(i) = meany+(covXY-varX)*(X(i)-meanx);
    ylSE(i) = (Y(i)-ylEstimator(i))^2;
end
% mean square error for linear estimator
LMMSE = sum(ylSE)/N;

%haven't calculated the theoretical value yet
T = table(MMSE, LMMSE);
display(T);