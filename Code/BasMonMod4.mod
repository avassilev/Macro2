/*
This file implements a version of 
the basic classical monetary model 
in Gali (2008), Chapter 2
augmented with a simple 
interest rate rule
*/

var a, c, p, m, i_a, n, i, Pi, r, omega, y; 
varexo eps_a
       eps_i;  

parameters alppha betta rho sigma varphi rho_a rho_i phi_pi eta;

alppha = 0.33; 
betta = 0.99;  
rho = -ln(betta);  
sigma = 0.5;    
varphi = 1;  
rho_a = 0.9;
rho_i = 0.9;
phi_pi = 0.5;   // Try changing this to a value between 0 and 1
eta = 4;



model(linear);
omega = sigma*c + varphi*n;               // Eq. (1)
c = c(+1) - 1/sigma*(r - rho);            // Eq. (2)
m - p = y - eta*i;                        // Eq. (3)
y = a + (1-alppha)*n;                     // Eq. (4)
omega = a - alppha*n + ln(1-alppha);      // Eq. (5)
c = y;                                    // Eq. (6)
a = rho_a*a(-1) + eps_a;                  // AR(1) process for log-technology
i = r + Pi(+1);
Pi = p - p(-1);

i = rho + phi_pi*Pi + i_a;                // Simple interest rate rule (modified with a shock)
i_a = rho_i*i_a(-1) + eps_i;
end;

shocks;
// var eps_a; stderr 1;
var eps_i; stderr 1;
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
p=5;
m = y + p - eta*i; 
i_a = 0;
end;

resid(1);
steady;
check;

