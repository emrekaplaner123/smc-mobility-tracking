function sis = run_sis(model, Y, pos_vec, N, hist_times)

Phi = model.Phi;
Psi_z = model.Psi_z;
Psi_w = model.Psi_w;
P = model.P;
mu0 = model.mu0;
Sigma0 = model.Sigma0;
Z_vals = model.Z_vals;
sigma = model.sigma;
v = model.v;
eta = model.eta;
varsigma = model.varsigma;

m = size(Y,2) - 1;

hist_times = hist_times(hist_times <= m);

cdfP = cumsum(P,2);
L = chol(Sigma0, 'lower');

X_part = randn(N,6) * L' + mu0';
z_idx = randi(5,N,1);

tau1 = zeros(m+1,1);
tau2 = zeros(m+1,1);
ESS = zeros(m+1,1);
histW = zeros(N,length(hist_times));

dx = X_part(:,1) - pos_vec(1,:);
dy = X_part(:,4) - pos_vec(2,:);
dist = sqrt(dx.^2 + dy.^2);
dist = max(dist,1e-12);

mu_y = v - 10*eta*log10(dist);
res = Y(:,1)' - mu_y;

logw = -0.5*sum((res./varsigma).^2,2) - 6*log(sqrt(2*pi)*varsigma);

a = max(logw);
w = exp(logw - a);
w = w / sum(w);

tau1(1) = sum(w .* X_part(:,1));
tau2(1) = sum(w .* X_part(:,4));
ESS(1) = 1 / sum(w.^2);

for j = 1:length(hist_times)
    if hist_times(j) == 0
        histW(:,j) = w;
    end
end

for n = 1:m
    Z = Z_vals(:, z_idx)';
    W = sigma * randn(N,2);

    X_part = X_part * Phi' + Z * Psi_z' + W * Psi_w';

    U = rand(N,1);
    z_idx = 1 + sum(bsxfun(@gt, U, cdfP(z_idx,:)), 2);

    dx = X_part(:,1) - pos_vec(1,:);
    dy = X_part(:,4) - pos_vec(2,:);
    dist = sqrt(dx.^2 + dy.^2);
    dist = max(dist,1e-12);

    mu_y = v - 10*eta*log10(dist);
    res = Y(:,n+1)' - mu_y;

    logw = logw + (-0.5*sum((res./varsigma).^2,2) ...
        - 6*log(sqrt(2*pi)*varsigma));

    a = max(logw);
    w = exp(logw - a);
    w = w / sum(w);

    tau1(n+1) = sum(w .* X_part(:,1));
    tau2(n+1) = sum(w .* X_part(:,4));
    ESS(n+1) = 1 / sum(w.^2);

    for j = 1:length(hist_times)
        if hist_times(j) == n
            histW(:,j) = w;
        end
    end
end

sis.tau1 = tau1;
sis.tau2 = tau2;
sis.ESS = ESS;
sis.histW = histW;
sis.hist_times = hist_times;
sis.m = m;

end