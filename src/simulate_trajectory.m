function X = simulate_trajectory(model, m)

Phi = model.Phi;
Psi_z = model.Psi_z;
Psi_w = model.Psi_w;
P = model.P;
mu0 = model.mu0;
Sigma0 = model.Sigma0;
Z_vals = model.Z_vals;
sigma = model.sigma;

L = chol(Sigma0, 'lower');
X_0 = mu0 + L*randn(6,1);

X = zeros(6,m+1);
X(:,1) = X_0;

z_index = randi(5);

for n = 1:m
    Z = Z_vals(:, z_index);
    W = sigma * randn(2,1);

    X(:,n+1) = Phi * X(:,n) + Psi_z * Z + Psi_w * W;

    u = rand;
    cdf_row = cumsum(P(z_index,:));
    z_index = find(u <= cdf_row, 1, 'first');
end

end