/*
This file implements a version of 
the basic classical monetary model 
in Gali (2008), Chapter 2
for the case where the central bank
sets the money supply
*/

var a, c, p, m, n, i, Pi, r, omega, y, u; 
varexo eps_a
       eps_m;  

parameters alppha betta rho sigma varphi rho_a rho_m eta;

alppha = 0.33; 
betta = 0.99;  
rho = -ln(betta);  
sigma = 0.5;    
varphi = 1;  
rho_a = 0.9;
rho_m = 0.9;
eta = 4;



model(linear);
omega = sigma*c + varphi*n;               // Eq. (1)
c = c(+1) - 1/sigma*(r - rho);            // Eq. (2)
i = 1/eta*(y - (m - p));                  // Eq. (3)
y = a + (1-alppha)*n;                     // Eq. (4)
omega = a - alppha*n + ln(1-alppha);      // Eq. (5)
c = y;                                    // Eq. (6)
a = rho_a*a(-1) - eps_a;                  // AR(1) process for log-technology

Pi = p - p(-1);
p = (eta/(1+eta))*p(+1) + (1/(1+eta))*m + u; // Eq. (15)
m = rho_m*m(-1) + eps_m;
u = (eta*r-y)/(1+eta);
end;

shocks;
// var eps_a; stderr 1;
var eps_m; stderr 1;
end;

initval;
a = 0;
r = rho;
n = ln(1-alppha)/(sigma*(1-alppha) + varphi + alppha);
y = (1-alppha)*n;
c = y;
omega = sigma*c + varphi*n;
Pi=0.0;
p=5;  // Try changing this value to see if it will affect the steady state
m = 0; 
i = 1/eta*(y - (m - p));
u = (eta*r-y)/(1+eta);
end;

resid(1);
steady;
check;

stoch_simul(irf=20, noprint);
