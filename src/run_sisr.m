function sisr = run_sisr(model, Y, pos_vec, N)

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

cdfP = cumsum(P,2);
L = chol(Sigma0, 'lower');

X_part = randn(N,6) * L' + mu0';
z_idx = randi(5,N,1);

tau1 = zeros(m+1,1);
tau2 = zeros(m+1,1);
cmd_prob = zeros(m+1,5);
cmd_mode = zeros(m+1,1);

dx = X_part(:,1) - pos_vec(1,:);
dy = X_part(:,4) - pos_vec(2,:);
dist = sqrt(dx.^2 + dy.^2);
dist = max(dist,1e-12);

mu_y = v - 10*eta*log10(dist);
res = Y(:,1)' - mu_y;

logg = -0.5*sum((res./varsigma).^2,2) - 6*log(sqrt(2*pi)*varsigma);

a = max(logg);
g = exp(logg - a);
w = g / sum(g);

tau1(1) = sum(w .* X_part(:,1));
tau2(1) = sum(w .* X_part(:,4));

for j = 1:5
    cmd_prob(1,j) = sum(w(z_idx == j));
end

[~, cmd_mode(1)] = max(cmd_prob(1,:));

ind = systematic_resampling(w);
X_part = X_part(ind,:);
z_idx = z_idx(ind);

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

    logg = -0.5*sum((res./varsigma).^2,2) - 6*log(sqrt(2*pi)*varsigma);

    a = max(logg);
    g = exp(logg - a);
    w = g / sum(g);

    tau1(n+1) = sum(w .* X_part(:,1));
    tau2(n+1) = sum(w .* X_part(:,4));

    for j = 1:5
        cmd_prob(n+1,j) = sum(w(z_idx == j));
    end

    [~, cmd_mode(n+1)] = max(cmd_prob(n+1,:));

    ind = systematic_resampling(w);
    X_part = X_part(ind,:);
    z_idx = z_idx(ind);
end

sisr.tau1 = tau1;
sisr.tau2 = tau2;
sisr.cmd_prob = cmd_prob;
sisr.cmd_mode = cmd_mode;
sisr.m = m;

end