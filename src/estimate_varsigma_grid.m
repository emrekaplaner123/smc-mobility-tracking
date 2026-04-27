function calib = estimate_varsigma_grid(model, Y_unknown, pos_vec, varsigma_grid, N)

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

m_unknown = size(Y_unknown,2) - 1;

loglik_grid = zeros(length(varsigma_grid),1);

cdfP = cumsum(P,2);
L = chol(Sigma0, 'lower');

for s = 1:length(varsigma_grid)

    varsigma_try = varsigma_grid(s);

    X_part = randn(N,6) * L' + mu0';
    z_idx = randi(5,N,1);

    dx = X_part(:,1) - pos_vec(1,:);
    dy = X_part(:,4) - pos_vec(2,:);
    dist = sqrt(dx.^2 + dy.^2);
    dist = max(dist,1e-12);

    mu_y = v - 10*eta*log10(dist);
    res = Y_unknown(:,1)' - mu_y;

    logg = -0.5*sum((res./varsigma_try).^2,2) ...
        - 6*log(sqrt(2*pi)*varsigma_try);

    a = max(logg);
    g = exp(logg - a);
    w = g / sum(g);

    loglik = a + log(mean(g));

    ind = systematic_resampling(w);
    X_part = X_part(ind,:);
    z_idx = z_idx(ind);

    for n = 1:m_unknown
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
        res = Y_unknown(:,n+1)' - mu_y;

        logg = -0.5*sum((res./varsigma_try).^2,2) ...
            - 6*log(sqrt(2*pi)*varsigma_try);

        a = max(logg);
        g = exp(logg - a);
        w = g / sum(g);

        loglik = loglik + a + log(mean(g));

        ind = systematic_resampling(w);
        X_part = X_part(ind,:);
        z_idx = z_idx(ind);
    end

    loglik_grid(s) = loglik / m_unknown;
end

[~, best_idx] = max(loglik_grid);
varsigma_hat = varsigma_grid(best_idx);

X_part = randn(N,6) * L' + mu0';
z_idx = randi(5,N,1);

tau1_hat = zeros(m_unknown+1,1);
tau2_hat = zeros(m_unknown+1,1);

dx = X_part(:,1) - pos_vec(1,:);
dy = X_part(:,4) - pos_vec(2,:);
dist = sqrt(dx.^2 + dy.^2);
dist = max(dist,1e-12);

mu_y = v - 10*eta*log10(dist);
res = Y_unknown(:,1)' - mu_y;

logg = -0.5*sum((res./varsigma_hat).^2,2) ...
    - 6*log(sqrt(2*pi)*varsigma_hat);

a = max(logg);
g = exp(logg - a);
w = g / sum(g);

tau1_hat(1) = sum(w .* X_part(:,1));
tau2_hat(1) = sum(w .* X_part(:,4));

ind = systematic_resampling(w);
X_part = X_part(ind,:);
z_idx = z_idx(ind);

for n = 1:m_unknown
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
    res = Y_unknown(:,n+1)' - mu_y;

    logg = -0.5*sum((res./varsigma_hat).^2,2) ...
        - 6*log(sqrt(2*pi)*varsigma_hat);

    a = max(logg);
    g = exp(logg - a);
    w = g / sum(g);

    tau1_hat(n+1) = sum(w .* X_part(:,1));
    tau2_hat(n+1) = sum(w .* X_part(:,4));

    ind = systematic_resampling(w);
    X_part = X_part(ind,:);
    z_idx = z_idx(ind);
end

calib.varsigma_grid = varsigma_grid;
calib.loglik_grid = loglik_grid;
calib.varsigma_hat = varsigma_hat;
calib.tau1_hat = tau1_hat;
calib.tau2_hat = tau2_hat;
calib.m = m_unknown;

end