/*
This file implements a version of 
the New Keynesian model 
in Gali (2008), Chapter 3,
with a simple interest rate rule
*/

var a, Pi, y_tilde, y, y_n, n, i, r, r_n, l, l_gr_ann, nu, Pi_ann, i_ann, r_ann; 
varexo eps_a
       eps_nu; 

parameters alppha betta rho sigma epsilon varphi mu rho_a rho_nu phi_pi phi_y eta theta Theta lambda kappa psi_n_ya vartheta_n_y;

alppha = 0.33; 
betta = 0.99;  
rho = -ln(betta);  
sigma = 1;    
epsilon = 6;
varphi = 1;  
mu = ln(epsilon/(epsilon-1));
rho_a = 0.9;
rho_nu = 0.5;
phi_pi = 1.5; 
phi_y = 0.5/4;
eta = 4;
theta = 2/3;
Theta = (1-alppha)/(1-alppha+alppha*epsilon);
lambda = (((1-theta)*(1-betta*theta))/theta)*Theta;
kappa = lambda*(sigma + (varphi+alppha)/(1-alppha));
psi_n_ya = (1+varphi)/(sigma*(1-alppha)+varphi+alppha);
vartheta_n_y = -((1-alppha)*(mu-ln(1-alppha)))/(sigma*(1-alppha)+varphi+alppha);


model(linear);
Pi = betta*Pi(+1) + kappa*y_tilde; // Eq. (20) 
y_tilde = -1/sigma*(i - Pi(+1) - r_n) + y_tilde(+1); // Eq. (21)
r_n = rho + sigma*psi_n_ya*(a(+1)-a); // Eq. (22)
y_n = psi_n_ya*a + vartheta_n_y; // Eq. (18)
y_tilde = y - y_n; // Output gap definition
l = y - eta*i; // Eq. (4) with l := m - p
l_gr_ann = 4*(l - l(-1)); // Annualized growth rate
y = a + (1-alppha)*n; // Eq. (12)
i = r + Pi(+1); // Fisher equation
a = rho_a*a(-1) + eps_a; // AR(1) process for log-technology
nu = rho_nu*nu(-1) + eps_nu; // AR(1) process for monetary policy shock
i = rho + phi_pi*Pi + phi_y*y_tilde + nu; // Eq. (24)
Pi_ann = 4*Pi; // Annualized inflation rate
i_ann = 4*i; // Annualized nominal interest rate
r_ann = 4*r; // Annualized real interest rate
end;

initval;
a = 0;
nu = 0;
r_n = rho;
r = rho;
y_n = vartheta_n_y;
y_tilde = 0;
y = y_n;
n = y/(1-alppha);
i=rho;
Pi=0.0;
l = y - eta*i; 
l_gr_ann = 0;
Pi_ann = 4*Pi;
i_ann = 4*i;
r_ann = 4*r;
end;

resid;
steady;
check;

shocks;
var eps_nu; stderr 0.25;
end;

stoch_simul(order = 1,irf=12) y_tilde Pi_ann i_ann r_ann l_gr_ann nu;

