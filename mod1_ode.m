function xdot = mod1_ode(~,y,q)
% The initial populations are set
S_b = y(1);
E_b = y(2);
I_b = y(3);
D_n = y(4)
S_n = y(5);
E_n = y(6);
I_n = y(7);

% We set the base params
beta_b = q(1);
alpha = q(2);
gamma_b = q(3);
mu_b = q(4);

beta_n = q(5);
gamma_n = q(6);
mu_n = q(7);

% Additional params
% None for the base model

% Bullfrog odes
S_b_dot = -beta_b * S_b * (E_b + I_b) - beta_b * S_b * (E_n + I_n) + alpha * I_b;
E_b_dot = beta_b * S_b * (E_b + I_b) + beta_b * S_b * (E_n + I_n) - gamma_b * E_b;
I_b_dot = gamma_b * E_b  - mu_b * I_b - alpha * I_b;

D_dot = mu_n * I_n;

% Native frog odes
S_n_dot =  - beta_n * S_n * (E_b + I_b) - beta_n * S_n * (E_n + I_n);
E_n_dot = beta_n * S_n * (E_b + I_b) + beta_n * S_n * (E_n + I_n) - gamma_n * E_n;
I_n_dot = gamma_n * E_n - mu_n * I_n;

xdot = [S_b_dot; E_b_dot; I_b_dot; D_dot; S_n_dot; E_n_dot; I_n_dot];

end