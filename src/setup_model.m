function model = setup_model()

dt = 0.5;
alpha = 0.6;
sigma = 0.5;

Phi_tilde = [1 dt dt^2/2;
             0 1 dt;
             0 0 alpha];

Psi_z_tilde = [dt^2/2;
               dt;
               0];

Psi_w_tilde = [dt^2/2;
               dt;
               1];

Phi = [Phi_tilde zeros(3,3);
       zeros(3,3) Phi_tilde];

Psi_z = [Psi_z_tilde zeros(3,1);
         zeros(3,1) Psi_z_tilde];

Psi_w = [Psi_w_tilde zeros(3,1);
         zeros(3,1) Psi_w_tilde];

P = (1/20) * [16 1 1 1 1;
              1 16 1 1 1;
              1 1 16 1 1;
              1 1 1 16 1;
              1 1 1 1 16];

mu0 = zeros(6,1);
Sigma0 = diag([500, 5, 5, 200, 5, 5]);

Z_vals = [0    3.5   0    0   -3.5;
          0    0     3.5 -3.5  0];

model.dt = dt;
model.alpha = alpha;
model.sigma = sigma;
model.Phi = Phi;
model.Psi_z = Psi_z;
model.Psi_w = Psi_w;
model.P = P;
model.mu0 = mu0;
model.Sigma0 = Sigma0;
model.Z_vals = Z_vals;
model.v = 90;
model.eta = 3;
model.varsigma = 1.5;

end