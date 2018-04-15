%% --- What is this script? ---
% This script solves the last question of the 1st exercise: Noise
% canceling.

%% --- Load data ---
load sounds.mat

%% --- Samples ---
N = size(d, 1);
n = 1:N;

%% --- Filter parameters ---
% Number of w params
K = 500;
% steepest descent scalar parameter
mu = 1.5;

%% --- Filter parameters and output ---
[R, p] = computePR(u, d, K);
[y, w, wt] = gradientDescent(u, R, p, mu);

%% --- Optimal filter parameters w (wo- Wiener solution) ---
wo = R \ p;

%% --- Best possible filter output (offline) ---
y_best = zeros(N, 1);
for i = K:N
    y_best(i) = sum(wo'*u(i:-1:i-K+1));
end

%% --- Output e (with the adaptive filter) ---
e = d - y;
