/*
This file implements a version of 
the basic classical monetary model 
in Gali (2008), Chapter 2
augmented with an exogenous
process for the nominal interest rate
*/

var a, c, p, m, n, i, Pi, r, omega, y; 
varexo eps_a
       eps_i;   

parameters alppha betta rho sigma varphi rho_a rho_i eta;

alppha = 0.33; 
betta = 0.99;  
rho = -ln(betta);  
sigma = 0.5;    
varphi = 1;  
rho_a = 0.9;
rho_i = 0.9;
eta = 4;



model(linear);
omega = sigma*c + varphi*n;               // Eq. (1)
c = c(+1) - 1/sigma*(r - rho);            // Eq. (2)
m - p = y - eta*i;                        // Eq. (3)
y = a + (1-alppha)*n;                     // Eq. (4)
omega = a - alppha*n + ln(1-alppha);      // Eq. (5)
c = y;                                    // Eq. (6)
a = rho_a*a(-1) + eps_a;                  // AR(1) process for log-technology
i = r + Pi(+1);                           // Eq. (11)
Pi = p - p(-1);

i = rho*(1-rho_i) + rho_i*i(-1) + eps_i;  // Exogenous process for the nominal interest rate
end;

initval;
a = 0;
r = rho;
n = ln(1-alppha)/(sigma*(1-alppha) + varphi + alppha);
y = (1-alppha)*n;
c = y;
omega = sigma*c + varphi*n;
i=rho+Pi;
Pi=0.0;
p=3; // Try changing this value to see if it will affect the steady state
m = y + p - eta*i; 
end;

resid(1);
steady;
check;

