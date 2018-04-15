%% --- What is this script? ---
% This script solves the first 3 questions of the 1st exercise: Noise
% canceling.

%% --- Samples ---
N = 50000;
n = 1:N;

%% --- Filter parameters ---
% Number of w params
K = 2;
% steepest descent parameters
mu = 1.5;
error_threshold = 1e-12;

%% --- x signal (clean signal) ---
x = cos(n*pi).*sin(pi/25*n+pi/3); x = x(:);

%% --- white noise ---
vr = 0.19;
v = sqrt(vr)*randn(N,1); v = v-mean(v);

%% --- noised signal ---
d = x+v;

%% --- Measured noise (noise with noise :P) ---
u = zeros(N, 1);
for i = 2:N
    u(i) = -0.78*u(i-1)+v(i);
end

%% --- Filter parameters and output ---
[R, p] = computePR(u, d, K);
[y, w, wt] = gradientDescent(u, R, p, mu);
%[y, w, wt] = gradientDescent(u, R, p, mu, error_threshold);
%% --- Optimal filter parameters w (wo- Wiener solution) ---
wo = R \ p;

%% --- Best possible filter output (offline) ---
y_best = zeros(N, 1);
for i = K:N
    y_best(i) = sum(w'*flip(u(i-K+1:i)));
end

%% --- Output e (with the adaptive filter) ---
e = d - y;

%% --- Output e_best (with the optimal filter) ---
e_best = d - y_best;

%% --- Max mu ---
mu_max = 2/max(eigs(R));

%% --- Contour plot of the error surface ---
L = 50;
WW = linspace(-2, 2, L);
J = zeros([L,L]);
for i = 1:L
    for j = 1:L
        WT = [WW(i); WW(j)];
        J(j, i) = var(d) - 2*p'*WT + WT'*R*WT;
    end
end

min_J = min(J(:));
max_J = max(J(:));

levels = linspace(min_J,max_J,20);

figure('name', 'J surface and adaptation process');
contour(WW, WW, J, levels/10); axis square
hold on
plot(wt(1,1), wt(2,1), '*');
plot(wt(1, :), wt(2, :), 'r-');
plot(wo(1), wo(2), 'ob')
hold off
colorbar
xlabel('w(2)');
ylabel('w(1)');
title('Error Surface and Adaptation process');

%% --- Jw criterion values ---
Jwt = zeros(N, 1);
sigma_d = var(d);
for i = 1:N
    Jwt(i) = sigma_d - 2*p'*wt(:, i) + wt(:, i)'*R*wt(:, i);
end

figure('name', 'Jw value after each step')
semilogx(n, Jwt);
title('Jw value after each step')
xlabel('step')
ylabel('J_w value')
