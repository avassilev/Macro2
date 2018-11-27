function [k,kplus] = OLG_CobbDouglas()
alppha = 0.33;
g = 0.1;
n = 0.05;
rho = 0.1;
kmax = 0.3;
step = 0.001;
    function kp = kprime(k)
    kp=(1-alppha)/((1+g)*(1+n)*(2+rho))*k.^alppha;
    end
k = 0:step:kmax;
kplus = kprime(k);
plot(k,k,'Color','black','LineWidth',2);
hold on;
plot(k,kplus,'Color','blue','LineWidth',2);
rho = 1;
kplus = kprime(k);
plot(k,kplus,'Color','green','LineWidth',2);
rho = 2;
kplus = kprime(k);
plot(k,kplus,'Color','red','LineWidth',2);
xlabel('$k_t$','Interpreter','latex','FontSize',30)
ylabel('$k_{t+1}$','Interpreter','latex','FontSize',30)
end