% This script should be run once the final values have been
% determined for all parameters. It will generate all of 
% our figures to include in our analysis. 

% Set upper_range to the time span you want
upper_range = 60;
tspan = 0:1:upper_range;

% --------- INITIAL POPULATIONS SET UP ---------

% The bullfrog populations
S_b0 = 8;
E_b0 = 1;
I_b0 = 1;

% The native frog populations
S_n0 = 100;
E_n0 = 0;
I_n0 = 0;

D_n0 = 0;

% Total frogs
total_frogs = S_b0 + E_b0 + I_b0 + S_n0 + E_n0 + I_n0;

x01 = [S_b0; E_b0; I_b0; D_n0; S_n0; E_n0; I_n0];
x02 = [S_b0; E_b0; I_b0; D_n0; S_n0; E_n0; I_n0];
x03 = [S_b0; E_b0; I_b0; D_n0; S_n0; E_n0; I_n0];
x04 = [S_b0; E_b0; I_b0; D_n0; S_n0; E_n0; I_n0];
x05 = [S_b0; E_b0; I_b0; D_n0; S_n0; E_n0; I_n0];

% --------- PARAM SET UP ---------

% --------- Constant across models ---------
% Set up bullfrog params
% Beta is the rate of transition for susceptible bullfrogs to exposed
beta_b = 0.005;

% Alpha is rate at which bullfrogs recovery from chytrid infection
alpha = 0.75;

% Rate of transition from exposed to infected for bullfrogs
gamma_b = 10/100;

% Death rate of infected bullfrogs from chytrid
mu_b = 1/200;


% Set up native frog params
% Beta is the rate of transition for susceptible frogs to exposed
beta_n = 0.01;

% Rate of transition from exposed to infected for native frogs
gamma_n = 70/100;

% Death rate of infected native frogs from chytrid
mu_n = 80/100; 

% --------- Model dependent ---------
% For removal model
% Rate of bull frog removal
sigma = 0.30;

% For washing models
% Rate of washing of natives
omega = 0.35;

% Make parameter lists
parms1 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n];
parms2 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n, sigma];
parms3 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n, omega];
parms4 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n, sigma, omega];


% --------- MODEL GENERATION ---------

% Run all of our odes
[~,X1]=ode23(@(t,x) mod1_ode(t,x,parms1),tspan,x01);
[~,X2]=ode23(@(t,x) mod2_ode(t,x,parms2),tspan,x02);
[~,X3]=ode23(@(t,x) mod3_ode(t,x,parms3),tspan,x03);
[~,X4]=ode23(@(t,x) mod4_ode(t,x,parms4),tspan,x04);

% -----------------------------------------------------
% CHART 1: Native frogs left living chart 
% Creating parameter of dead natives in the presence of bullfrogs
n_dead = [X1(1:(upper_range+1), 4) X2(1:(upper_range+1), 4) X3(1:(upper_range+1), 4) X4(1:(upper_range+1), 4)];

%disp(mat2str(size(n_dead)));

% Create matrix with all 100 entries
A = zeros(upper_range+1,4) + 100;

%disp(mat2str(size(A)));

% Get the complement of the dead and total frogs to get
% the frog population still alive
n_alive = A - n_dead;

figure(1);
set(gcf,'color','white')
plot(tspan,n_alive)

ylim([0 total_frogs])
leg = legend('Base model', 'Remove Bullfrogs', 'Wash Natives','Native Washing and Bullfrog Removal');
leg.Location = 'northeastoutside';
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Living Native Population')
title('Living Natives over Time')
% -----------------------------------------------------
% CHART 2: Bullfrogs left alive

b_alive = [sum(X1(:,[2 3]),2) sum(X2(:,[2 3]),2) sum(X3(:,[2 3]),2) sum(X4(:,[2 3]),2)];

figure(2);
set(gcf,'color','white')
plot(tspan,b_alive)

ylim([0 12])
leg = legend('Base model', 'Remove Bullfrogs', 'Wash Natives','Native Washing and Bullfrog Removal');
leg.Location = 'northeastoutside';
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Living Bullfrog Population')

title('Living Bullfrogs over Time')
% -----------------------------------------------------
% CHART 3: Optimal sigma (living frog population with changing sigma)

%Simple for-loop analysis, do it with changing parameter.

%Omega constant at 0.1
omega = 0.134;

%sigma = 0.184;

%omega = 0.051, sigma = 0.096
%omega = 0.134, sigma = 0.184

max_n_alive_model_four = 0;

sigma_vals = 0:0.001:0.25;

alive_vals = [];

for sigma_val = 0:0.001:0.25
    
sigma = sigma_val;

parms1 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n];
parms2 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n, sigma];
parms3 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n, omega];
parms4 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n, sigma, omega];

[~,X1]=ode23(@(t,x) mod1_ode(t,x,parms1),tspan,x01);
[~,X2]=ode23(@(t,x) mod2_ode(t,x,parms2),tspan,x02);
[~,X3]=ode23(@(t,x) mod3_ode(t,x,parms3),tspan,x03);
[~,X4]=ode23(@(t,x) mod4_ode(t,x,parms4),tspan,x04);

n_dead = [X1(1:(upper_range+1), 4) X2(1:(upper_range+1), 4) X3(1:(upper_range+1), 4) X4(1:(upper_range+1), 4)];

A = zeros(upper_range+1,4) + 100;

n_alive = A - n_dead;

n_alive_model_one = n_alive(61, 1);
n_alive_model_two = n_alive(61, 2);
n_alive_model_three = n_alive(61, 3);
n_alive_model_four = n_alive(61, 4);

alive_vals(end + 1) = n_alive_model_four;

%break

%end

end

figure(3);
set(gcf,'color','white')

display(mat2str(size(sigma_vals)))
display(mat2str(size(alive_vals)))

plot(sigma_vals,alive_vals)

ylim([0 total_frogs])
leg = legend('Native Washing and Bullfrog Removal');
leg.Location = 'northeastoutside';
set(gca,'fontsize',16)
xlabel('Sigma')
ylabel('Living Native Population')
title('Living Natives over different Sigma values')

% -----------------------------------------------------
% CHART 4: Optimal omega (living frog population with changing omega)

%Simple for-loop analysis, do it with changing parameter.

%Main logic.
sigma = 0.184;
%omega = 0.184;

omega_vals = 0:0.001:0.25;

alive_vals = [];

%max_n_alive_model_four = 0;

for omega_val = 0:0.001:0.25
    
omega = omega_val;

parms1 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n];
parms2 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n, sigma];
parms3 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n, omega];
parms4 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n, sigma, omega];

[~,X1]=ode23(@(t,x) mod1_ode(t,x,parms1),tspan,x01);
[~,X2]=ode23(@(t,x) mod2_ode(t,x,parms2),tspan,x02);
[~,X3]=ode23(@(t,x) mod3_ode(t,x,parms3),tspan,x03);
[~,X4]=ode23(@(t,x) mod4_ode(t,x,parms4),tspan,x04);

n_dead = [X1(1:(upper_range+1), 4) X2(1:(upper_range+1), 4) X3(1:(upper_range+1), 4) X4(1:(upper_range+1), 4)];

A = zeros(upper_range+1,4) + 100;

n_alive = A - n_dead;

n_alive_model_one = n_alive(61, 1);
n_alive_model_two = n_alive(61, 2);
n_alive_model_three = n_alive(61, 3);
n_alive_model_four = n_alive(61, 4);

alive_vals(end+1) = n_alive_model_four;

end

figure(4);
set(gcf,'color','white')

display(mat2str(size(omega_vals)))
display(mat2str(size(alive_vals)))

plot(omega_vals,alive_vals)

ylim([0 total_frogs])
leg = legend('Native Washing and Bullfrog Removal');
leg.Location = 'northeastoutside';
set(gca,'fontsize',16)
xlabel('Omega')
ylabel('Living Native Population')
title('Living Natives over different Omega values')


