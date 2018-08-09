% Vectorized Monte Carlo implementation of evaluating an asian option 
% with the Milstein sheme
% Computaional Finance, University of St. Gallen, 2017

tic
M           = 100000;  %Number of paths          
S0          = 100;     %Initital Asset price
sigma20     = 0.25;    %Initial vol  
T           = 1;       %Time to maturity
rho         = 0.05;    %Interest rates
kappa       = 2;       %Rate of mean reversion (vol)
theta       = 0.4;     %Long term mean (vol)
nu          = 0.2;     %volatility parameter of the variance process 0.2
p           = 0.2;     %correlation asset price and vol

N           = 100;     %number of steps
h           = T/N;     %Step size
 

randn('state',1)     % Reset Random Number Generator to 
  

%Initialize vectors and arrays
   
S = zeros(1,N);
V = zeros(1,N);   
S(1) = S0;
V(1) = sigma20;
   
avg_price  = zeros(1,M);
payoff     = zeros(1,M);
      
% 2 dimensional Wiener-Process

dW1 = randn(M,N+1)*sqrt(h);

dW2 = rho*dW1 + sqrt(1-rho^2)*randn(M,N+1)*sqrt(h);

% Intitialiszing asset price and sigma squared

S = S0*ones(M,N+1);

sigma2 = sigma20*ones(M,N+1);


% Slution of SDE system with the Milstein sheme

for i = 1:N
 
    sigma2(:,i+1) = sigma2(:,i) + kappa*(theta-sigma2(:,i))*h + nu*sqrt(sigma2(:,i)).*dW2(:,i) + (nu^2)/4*((dW2(:,i)).^2-h);
    
    S(:,i+1)      = S(:,i).*(1 + rho*h + sqrt(sigma2(:,i)).*dW1(:,i) + 0.5 .*sigma2(:,i).*((dW2(:,i)).^2-h));
    
end

payoff = max(0,S(:,N+1) - mean(S,2));

V = exp(-rho*T)*mean(payoff)
toc
