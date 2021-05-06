% Set upper_range to the time span you want
upper_range = 20;
tspan = 0:1:upper_range;

% --------- INITIAL POPULATIONS SET UP ---------

% The bullfrog populations
S_b0 = 10;
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
beta_b = 0.10;

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
gamma_n = 75/100;

% Death rate of infected native frogs from chytrid
mu_n = 80/100; 

% --------- Model dependent ---------
% For removal model
% Rate of bull frog removal
sigma = 0.20;

% For washing models
% Rate of washing of natives
omega = 0.20;

% Make parameter lists
parms1 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n];
parms2 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n, sigma];
parms3 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n, omega];
parms4 = [beta_b, alpha, gamma_b, mu_b, beta_n, gamma_n, mu_n, sigma, omega];


% --------- MODEL GENERATION ---------

% Run all of our odes
[~,X3]=ode23(@(t,x) mod3_ode(t,x,parms3),tspan,x03);
[~,X4]=ode23(@(t,x) mod4_ode(t,x,parms4),tspan,x04);

% Gen Model 1
[~,X1]=ode23(@(t,x) mod1_ode(t,x,parms1),tspan,x01);

figure(1);
set(gcf,'color','white')
plot(tspan,X1)
ylim([0 total_frogs])
leg = legend('S_b', 'E_b', 'I_b', 'D_n', 'S_n', 'E_n', 'I_n', 'Location','northeastoutside');
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Population')
title('Model 1 - Base Model - Solution')

% Gen Model 2
[~,X2]=ode23(@(t,x) mod2_ode(t,x,parms2),tspan,x02);

figure(2);
set(gcf,'color','white')
plot(tspan,X1)
ylim([0 total_frogs])
leg = legend('S_b', 'E_b', 'I_b', 'D_n', 'S_n', 'E_n', 'I_n', 'Location','northeastoutside');
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Population')
title('Model 2 - Remove Bullfrogs - Solution')

% Gen Model 3
[~,X3]=ode23(@(t,x) mod3_ode(t,x,parms3),tspan,x03);

figure(3);
set(gcf,'color','white')
plot(tspan,X1)
ylim([0 total_frogs])
leg = legend('S_b', 'E_b', 'I_b', 'D_n', 'S_n', 'E_n', 'I_n', 'Location','northeastoutside');
set(gca,'fontsize',16)
xlabel('Time (days)')
ylabel('Population')
title('Model 2 - Wash Natives - Solution')

