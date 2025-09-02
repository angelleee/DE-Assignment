clc; clear all; close all;
seed = 101;
rng(seed);

A0 = 2;
A = (2*A0) * rand - A0;
sigma2 = 0.1;
N = 10;
sigma2N = sigma2/N;

x = A + sigma2 * randn(N, 15);
x_ = mean(x, 1);

figure;
xlabel('A');
ylabel('P(A|x)');
title('Problem 2');
hold on;

for i = 1:15
    P_ax = @(a) (a >= -A0 & a <= A0) .* ( 1 / sqrt(2*pi*sigma2N) ) .* exp( (-0.5/sigma2N) * ((a-x_(i)).^2) );
    
    c = integral(P_ax, -A0, A0);
    P_ax_normalized = @(a) P_ax(a) / c;
    P_ax_mean = @(a) P_ax_normalized(a) .* a;
    A_mmse = integral(P_ax_mean, -A0, A0);
    
    fplot(P_ax_normalized, [-A0, A0]);
    plot(A_mmse, P_ax_normalized(A_mmse), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineWidth', 2);
end
hold off;

N = 100;
sigma2N = sigma2/N;

x = A + sigma2 * randn(N, 15);
x_ = mean(x, 1);

figure;
xlabel('A');
ylabel('P(A|x)');
title('Problem 3');
hold on;

for i = 1:15
    P_ax = @(a) (a >= -A0 & a <= A0) .* ( 1 / sqrt(2*pi*sigma2N) ) .* exp( (-0.5/sigma2N) * ((a-x_(i)).^2) );
    
    c = integral(P_ax, -A0, A0);
    P_ax_normalized = @(a) P_ax(a) / c;
    P_ax_mean = @(a) P_ax_normalized(a) .* a;
    A_mmse = integral(P_ax_mean, -A0, A0);
    
    fplot(P_ax_normalized, [-A0, A0]);
    plot(A_mmse, P_ax_normalized(A_mmse), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineWidth', 2);
end
hold off;

M = 1000;
A0 = 2;
sigma2 = 0.1;
N = 10;
sigma2N = sigma2/N;

A = (2*A0) * rand(1, M) - A0;

A_mmse = zeros(1,M);
for i = 1:M
    x = A(i) + sigma2 * randn(N, 1);
    x_ = mean(x);
    P_ax = @(a) (a >= -A0 & a <= A0) .* ( 1 / sqrt(2*pi*sigma2N) ) .* exp( (-0.5/sigma2N) * ((a-x_).^2) );
    c = integral(P_ax, -A0, A0);
    P_ax_normalized = @(a) P_ax(a) / c;
    P_ax_mean = @(a) P_ax_normalized(a) .* a;
    A_mmse(i) = integral(P_ax_mean, -A0, A0);
end
A_mmse_monte = mean(A_mmse);
BMSE = mean((A-A_mmse_monte).^2);
disp(['Problem(4) A_mmse using Monte Carlo: ', num2str(A_mmse_monte)]);
disp(['Problem(4) BMSE using Monte Carlo: ', num2str(BMSE)]);


N = 100;
sigma2N = sigma2/N;

A_mmse = zeros(1,M);
for i = 1:M
    x = A(i) + sigma2 * randn(N, 1);
    x_ = mean(x);
    P_ax = @(a) (a >= -A0 & a <= A0) .* ( 1 / sqrt(2*pi*sigma2N) ) .* exp( (-0.5/sigma2N) * ((a-x_).^2) );
    c = integral(P_ax, -A0, A0);
    P_ax_normalized = @(a) P_ax(a) / c;
    P_ax_mean = @(a) P_ax_normalized(a) .* a;
    A_mmse(i) = integral(P_ax_mean, -A0, A0);
end
A_mmse_monte = mean(A_mmse);
BMSE = mean((A-A_mmse_monte).^2);
disp(['Problem(5) A_mmse using Monte Carlo: ', num2str(A_mmse_monte)]);
disp(['Problem(5) BMSE using Monte Carlo: ', num2str(BMSE)]);

%%
seed = 101;
rng(seed);

A0 = 2;
A = (2*A0) * rand - A0;
sigma2 = 0.1;
N = 10;
sigma2N = sigma2/N;

x = A + sigma2 * randn(N, 15);
x_ = mean(x, 1);

figure;
xlabel('A');
ylabel('P(A|x)');
title('Problem 6-2');
hold on;

for i = 1:15
    P_ax = @(a) (a >= -A0 & a <= A0) .* ( 1 / sqrt(2*pi*sigma2N) ) .* exp( (-0.5/sigma2N) * ((a-x_(i)).^2) );
    
    c = integral(P_ax, -A0, A0);
    P_ax_normalized = @(a) P_ax(a) / c;
    CDF = @(a) integral(P_ax_normalized, -A0, a);
    median_A = fzero(@(a) CDF(a) - 0.5, 0);
    
    fplot(P_ax_normalized, [-A0, A0]);
    plot(median_A, P_ax_normalized(median_A), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineWidth', 2);
end
hold off;

N = 100;
sigma2N = sigma2/N;

x = A + sigma2 * randn(N, 15);
x_ = mean(x, 1);

figure;
xlabel('A');
ylabel('P(A|x)');
title('Problem 6-3');
hold on;

for i = 1:15
    P_ax = @(a) (a >= -A0 & a <= A0) .* ( 1 / sqrt(2*pi*sigma2N) ) .* exp( (-0.5/sigma2N) * ((a-x_(i)).^2) );
    
    c = integral(P_ax, -A0, A0);
    P_ax_normalized = @(a) P_ax(a) / c;
    CDF = @(a) integral(P_ax_normalized, -A0, a);
    median_A = fzero(@(a) CDF(a) - 0.5, 0);
    
    fplot(P_ax_normalized, [-A0, A0]);
    plot(median_A, P_ax_normalized(median_A), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineWidth', 2);
end
hold off;

M = 1000;
A0 = 2;
sigma2 = 0.1;
N = 10;
sigma2N = sigma2/N;

A = (2*A0) * rand(1, M) - A0;

median_A = zeros(1,M);
for i = 1:M
    x = A(i) + sigma2 * randn(N, 1);
    x_ = mean(x);
    P_ax = @(a) (a >= -A0 & a <= A0) .* ( 1 / sqrt(2*pi*sigma2N) ) .* exp( (-0.5/sigma2N) * ((a-x_).^2) );
    c = integral(P_ax, -A0, A0);
    P_ax_normalized = @(a) P_ax(a) / c;
    CDF = @(a) integral(P_ax_normalized, -A0, a);
    median_A(i) = fzero(@(a) CDF(a) - 0.5, 0);
end
A_hat = mean(median_A);
disp(['Problem(6-4) A_hat: ', num2str(A_hat)]);


N = 100;
sigma2N = sigma2/N;

median_A = zeros(1,M);
for i = 1:M
    x = A(i) + sigma2 * randn(N, 1);
    x_ = mean(x);
    P_ax = @(a) (a >= -A0 & a <= A0) .* ( 1 / sqrt(2*pi*sigma2N) ) .* exp( (-0.5/sigma2N) * ((a-x_).^2) );
    c = integral(P_ax, -A0, A0);
    P_ax_normalized = @(a) P_ax(a) / c;
    CDF = @(a) integral(P_ax_normalized, -A0, a);
    median_A(i) = fzero(@(a) CDF(a) - 0.5, 0);
end
A_hat = mean(median_A);
disp(['Problem(6-5) A_hat: ', num2str(A_hat)]);

%%
seed = 101;
rng(seed);

A0 = 2;
A = (2*A0) * rand - A0;
sigma2 = 0.1;
N = 10;
sigma2N = sigma2/N;

x = A + sigma2 * randn(N, 15);
x_ = mean(x, 1);

figure;
xlabel('A');
ylabel('P(A|x)');
title('Problem 7-2');
hold on;

for i = 1:15
    P_ax = @(a) (a >= -A0 & a <= A0) .* ( 1 / sqrt(2*pi*sigma2N) ) .* exp( (-0.5/sigma2N) * ((a-x_(i)).^2) );
    
    c = integral(P_ax, -A0, A0);
    P_ax_normalized = @(a) P_ax(a) / c;

    [max_a, max_val] = fminbnd(@(a) -P_ax_normalized(a), -A0, A0);
    
    fplot(P_ax_normalized, [-A0, A0]);
    plot(max_a, P_ax_normalized(max_a), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineWidth', 2);
end
hold off;

N = 100;
sigma2N = sigma2/N;

x = A + sigma2 * randn(N, 15);
x_ = mean(x, 1);

figure;
xlabel('A');
ylabel('P(A|x)');
title('Problem 7-3');
hold on;

for i = 1:15
    P_ax = @(a) (a >= -A0 & a <= A0) .* ( 1 / sqrt(2*pi*sigma2N) ) .* exp( (-0.5/sigma2N) * ((a-x_(i)).^2) );
    
    c = integral(P_ax, -A0, A0);
    P_ax_normalized = @(a) P_ax(a) / c;

    [max_a, max_val] = fminbnd(@(a) -P_ax_normalized(a), -A0, A0);
    
    fplot(P_ax_normalized, [-A0, A0]);
    plot(max_a, P_ax_normalized(max_a), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineWidth', 2);
end
hold off;

M = 1000;
A0 = 2;
sigma2 = 0.1;
N = 10;
sigma2N = sigma2/N;

A = (2*A0) * rand(1, M) - A0;

MAP_A = zeros(1,M);
for i = 1:M
    x = A(i) + sigma2 * randn(N, 1);
    x_ = mean(x);
    P_ax = @(a) (a >= -A0 & a <= A0) .* ( 1 / sqrt(2*pi*sigma2N) ) .* exp( (-0.5/sigma2N) * ((a-x_).^2) );
    c = integral(P_ax, -A0, A0);
    P_ax_normalized = @(a) P_ax(a) / c;
    [max_a, max_val] = fminbnd(@(a) -P_ax_normalized(a), -A0, A0);
    MAP_A(i) = max_a;
end
A_hat = mean(MAP_A);
disp(['Problem(7-4) A_hat: ', num2str(A_hat)]);


N = 100;
sigma2N = sigma2/N;

MAP_A = zeros(1,M);
for i = 1:M
    x = A(i) + sigma2 * randn(N, 1);
    x_ = mean(x);
    P_ax = @(a) (a >= -A0 & a <= A0) .* ( 1 / sqrt(2*pi*sigma2N) ) .* exp( (-0.5/sigma2N) * ((a-x_).^2) );
    c = integral(P_ax, -A0, A0);
    P_ax_normalized = @(a) P_ax(a) / c;
    [max_a, max_val] = fminbnd(@(a) -P_ax_normalized(a), -A0, A0);
    MAP_A(i) = max_a;
end
A_hat = mean(MAP_A);
disp(['Problem(7-5) A_hat: ', num2str(A_hat)]);
