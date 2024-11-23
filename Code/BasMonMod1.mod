/*
This file implements a version of 
the basic classical monetary model 
in Gali (2008), Chapter 2
*/

var a, c, n, r, omega, y; 
varexo eps_a;

parameters alppha betta rho sigma varphi rho_a;

alppha=0.33; 
betta=0.99;  
rho=-ln(betta);  
sigma=0.5;                                // Try changing this to 1.5 or 1.0 instead of 0.5
varphi=1;  
rho_a = 0.9;

model(linear);
omega = sigma*c + varphi*n;               // Eq. (1)
c = c(+1) - 1/sigma*(r - rho);            // Eq. (2)
y = a + (1-alppha)*n;                     // Eq. (4)
omega = a - alppha*n + ln(1-alppha);      // Eq. (5)
c = y;                                    // Eq. (6)
a = rho_a*a(-1) + eps_a;                  // AR(1) process for log-technology
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
end;

resid;
steady;
check;

stoch_simul(irf=30);