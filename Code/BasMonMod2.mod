/*
This file compares the results for 
two alternative implementations of 
the basic monetary model:
1) directly inputting the model 
   equations
2) using the explicit formulas 
   derived in Gali (2008), Chapter 2
*/

var a, c, n, n1, r, r1, omega, omega1, y, y1; 
varexo eps_a;

parameters alppha betta rho sigma varphi rho_a psi_na vartheta_n psi_ya vartheta_y psi_omegaa vartheta_omega;

alppha=0.33; 
betta=0.99;  
rho=-ln(betta);  
sigma=0.5;    
varphi=1;  
rho_a = 0.9;
psi_na = (1-sigma)/(sigma*(1-alppha) + varphi + alppha);
vartheta_n = ln(1-alppha)/(sigma*(1-alppha) + varphi + alppha);
psi_ya = (1+varphi)/(sigma*(1-alppha) + varphi + alppha);
vartheta_y = (1-alppha)*vartheta_n;
psi_omegaa = (sigma+varphi)/(sigma*(1-alppha) + varphi + alppha);
vartheta_omega = ((sigma*(1-alppha)+varphi)*ln(1-alppha))/(sigma*(1-alppha) + varphi + alppha);


model(linear);
omega = sigma*c + varphi*n;               // Eq. (1)
c = c(+1) - 1/sigma*(r - rho);            // Eq. (2)
y = a + (1-alppha)*n;                     // Eq. (4)
omega = a - alppha*n + ln(1-alppha);      // Eq. (5)
c = y;                                    // Eq. (6)
a = rho_a*a(-1) + eps_a;                  // AR(1) process for log-technology

n1 = psi_na*a + vartheta_n;               // Eq. (7)
y1 = psi_ya*a + vartheta_y;               // Eq. (8)
r1 = rho + sigma*psi_ya*(a(+1)-a);        // Eq. (9)
omega1 = psi_omegaa*a + vartheta_omega;   // Eq. (10)
end;

shocks;
var eps_a; stderr 1;
end;

initval;
a = 0;
r = rho;
n = ln(1-alppha)/(sigma*(1-alppha) + varphi + alppha);
y = (1-alppha)*n;
c = y;
omega = sigma*c + varphi*n;
n1 = psi_na*a + vartheta_n;
y1 = psi_ya*a + vartheta_y;
r1 = rho;
omega1 = psi_omegaa*a + vartheta_omega;
end;

resid(1);
steady;
check;

stoch_simul(irf=20, noprint, nograph);

diffs = [n_eps_a-n1_eps_a, y_eps_a-y1_eps_a, r_eps_a-r1_eps_a, omega_eps_a-omega1_eps_a]