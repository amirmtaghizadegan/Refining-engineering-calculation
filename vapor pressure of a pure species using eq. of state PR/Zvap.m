function Zv1=Zvap(beta,q,epsilon,sigma)
Zv1=1+beta-q*beta*(1-beta)/((1+epsilon*beta)*(1+sigma*beta));
err=1;
while err>eps
    Zv2=1+beta-q*beta*(Zv1-beta)/((Zv1+epsilon*beta)*(Zv1+sigma*beta));
    err=(Zv2-Zv1)/Zv2;
    Zv1=Zv2;
end
end