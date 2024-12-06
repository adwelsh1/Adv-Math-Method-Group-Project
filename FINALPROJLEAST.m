% SIR Model / Least Squares

clc; clear; close all;

N = 1000;        
S0 = 990;        
I0 = 10;          
R0 = 0;           
beta = 0.3;      
gamma = 0.1;     
h = 1;           
days = 30;       

% 30 days
t = 0:h:days;    
S = zeros(size(t)); I = zeros(size(t)); R = zeros(size(t));
S(1) = S0; I(1) = I0; R(1) = R0;


for k = 1:length(t)-1
    dS = -beta * S(k) * I(k) / N;
    dI = beta * S(k) * I(k) / N - gamma * I(k);
    dR = gamma * I(k);
    S(k+1) = S(k) + h * dS;
    I(k+1) = I(k) + h * dI;
    R(k+1) = R(k) + h * dR;
end

true_I = I;

%Linear least squares 

ln_I = log(true_I);     
k = 1:length(t);

% Least squares for 30 days 
X = [ones(length(t), 1), t(:)];  
y = ln_I(:);                    
theta = X \ y;                   
ln_I0_est = theta(1);            
k_est = theta(2);                
I0_est = exp(ln_I0_est);        
beta_est = (k_est + gamma) * N / S0;  

%10 days
t_10 = t(1:10);                  
ln_I_10 = log(true_I(1:10));    
X_10 = [ones(length(t_10), 1), t_10(:)];
y_10 = ln_I_10(:);
theta_10 = X_10 \ y_10;
ln_I0_est_10 = theta_10(1);
k_est_10 = theta_10(2);
I0_est_10 = exp(ln_I0_est_10);
beta_est_10 = (k_est_10 + gamma) * N / S0;


fprintf('--- Results using 30 days of data ---\n');
fprintf('Estimated I(0): %.2f\n', I0_est);
fprintf('Estimated beta: %.4f\n\n', beta_est);

fprintf('--- Results using 10 days of data ---\n'); 
fprintf('Estimated I(0): %.2f\n', I0_est_10);
fprintf('Estimated beta: %.4f\n\n', beta_est_10);


fprintf('Comparison of Estimates:\n');
fprintf('Difference in I(0): %.2f\n', abs(I0_est - I0_est_10));
fprintf('Difference in beta: %.4f\n', abs(beta_est - beta_est_10));


figure;
plot(t, true_I, 'b', 'LineWidth', 2); hold on;
plot(t, exp(theta(1) + theta(2) * t), 'r--', 'LineWidth', 2);  % Model fit for 30 days
plot(t_10, exp(theta_10(1) + theta_10(2) * t_10), 'g--', 'LineWidth', 2);  % Model fit for 10 days
xlabel('Time (days)');
ylabel('Infected Population');
legend('True I(t)', 'Fit (30 days)', 'Fit (10 days)');
title('SIR Model and Least Squares Fit');
grid on;
