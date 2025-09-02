A = 1;
M = 10000;
N = 2:1:100;
E_A = zeros(1,99);
var_A = zeros(1,99);

seed = 101;
rng(seed);

for j = N
    A_MLE = zeros(1,M);
    for i = 1:M
        x = normrnd(A,sqrt(A),1,j);
        sigma_x2 = sum(x.^2);
        A_MLE(i) = -1/2+sqrt(1/4+sigma_x2/j);
    end
    E_A(j-1) = sum(A_MLE)/M;
    A_MLE = A_MLE-E_A(j-1);
    var_A(j-1) = sum(A_MLE.^2)/M;
end

figure;
sgtitle('Problem (b)');

subplot(1,2,1);
plot(N,E_A);
xlabel('N');
ylabel('Monte Carlo mean');
title('Monte Carlo mean');

subplot(1,2,2);
plot(N,var_A);
xlabel('N');
ylabel('Monte Carlo variance');
title('Monte Carlo variance');

%%
A = 1;
N = [5,20,100];
M = 10000;

figure;
sgtitle('Problem (c)');
for n = 1:3
    j = N(n);
    A_MLE = zeros(1,M);
    for i = 1:M
        x = normrnd(A,sqrt(A),1,j);
        sigma_x2 = sum(x.^2);
        A_MLE(i) = -1/2+sqrt(1/4+sigma_x2/j);
    end

    x_Gaussian = 0:0.1:3;
    sigma_the = (A^2)/(j*(A+1/2));
    y_Gaussian = normpdf(x_Gaussian, A, sigma_the);

    subplot(1,3,n);
    histogram(A_MLE, 30, 'Normalization', 'pdf');
    hold on;
    plot(x_Gaussian, y_Gaussian);
    title(['N = ', num2str(j)]);
    legend('Histogram', 'Theoretical asymptotic pdf');
    xlabel('A');
    ylabel('probability');
end

%%
A = 1;
M = 10000;
N = 2;

seed = 101;
rng(seed);

MSE = zeros(1,1000);

while 1
    A_MLE = zeros(1,M);
    for i = 1:M
        x = normrnd(A,sqrt(A),1,N);
        sigma_x2 = sum(x.^2);
        A_MLE(i) = -1/2+sqrt(1/4+sigma_x2/N);
    end

    MSE(N-1) = mean((A_MLE-A).^2);

    if MSE(N-1) < 1e-3
        break
    end

    N = N+1;
end

figure;
plot(2:N,MSE(1:N-1));
hold on;
plot(N,MSE(N-1),'ro');
text(N-80,MSE(N-1)+0.02,['N_m_i_n(',num2str(N),',',num2str(MSE(N-1)),')']);
xlabel('N');
ylabel('MSE');
title('Problem (d)');

%%
A = 1;
N = [5,20,100];
M = 10000;
u = zeros(1,M);

g_double_prime = @(x) (-1/4).*(x+1/4).^(-3/2);
u0 = A+A^2;
error_bound_func = @(x) (abs(g_double_prime(x))/2) * (x-u0).^2;

figure;
sgtitle('Problem (e)');
for n = 1:3
    j = N(n);
    error = zeros(1,M);
    error_bound = zeros(1,M);
    for i = 1:M
        x = normrnd(A,sqrt(A),1,j);
        u(i) = sum(x.^2)/j;
        error_bound(i) = error_bound_func(u(i));
    end
    subplot(2,3,n);
    histogram(u, 50, 'Normalization', 'probability');
    title(['N = ', num2str(j), ', u']);
    xlabel('u');
    ylabel('u probability');

    subplot(2,3,n+3);
    histogram(error_bound, 50, 'Normalization', 'probability');
    title(['N = ', num2str(j), ', error bound']);
    xlabel('error bound');
    ylabel('error bound probability');
end

%%
A = 1;
N = [5,20,100];
M = 10000;
u = zeros(1,M);

g = @(x) (-1/2)+sqrt(1/4+x);
g_prime = @(x) (1/2).*(x+1/4).^(-1/2);
u0 = A+A^2;
Taylor_approx = @(x) g(u0)+g_prime(u0)*(x-u0);

figure;
sgtitle('Problem (f)');
for n = 1:3
    j = N(n);
    Taylor = zeros(1,M);
    for i = 1:M
        x = normrnd(A,sqrt(A),1,j);
        u(i) = sum(x.^2)/j;
        Taylor(i) = Taylor_approx(u(i));
    end

    x_Gaussian = 0:0.1:3;
    y_Gaussian = normpdf(x_Gaussian, mean(Taylor), var(Taylor));

    subplot(2,3,n);
    histogram(Taylor, 50, 'Normalization', 'pdf');
    title(['N = ', num2str(j)]);
    xlabel('A_M_L_E');
    ylabel('A_M_L_E pdf');

    subplot(2,3,n+3);
    plot(x_Gaussian, y_Gaussian);
    title(['N = ', num2str(j)]);
    xlabel('A_M_L_E');
    ylabel('Gaussian pdf');
end

%%
A = 1;
M = 10000;
N = 2;

g = @(x) (-1/2)+sqrt(1/4+x);
g_prime = @(x) (1/2).*(x+1/4).^(-1/2);
u0 = A+A^2;
Taylor_approx = @(x) g(u0)+g_prime(u0)*(x-u0);

seed = 101;
rng(seed);

MSE = zeros(1,1000);

while 1
    Taylor = zeros(1,M);
    u = zeros(1,M);

    for i = 1:M
        x = normrnd(A,sqrt(A),1,N);
        u(i) = sum(x.^2)/N;
        Taylor(i) = Taylor_approx(u(i));
    end
    
    MSE(N-1) = mean((Taylor-A).^2);

    if MSE(N-1) < 1e-3
        break
    end

    N = N+1;
end

figure;
plot(2:N,MSE(1:N-1));
hold on;
plot(N,MSE(N-1),'ro');
text(N-80,MSE(N-1)+0.02,['N_m_i_n(',num2str(N),',',num2str(MSE(N-1)),')']);
xlabel('N');
ylabel('MSE');
title('Problem (g)');