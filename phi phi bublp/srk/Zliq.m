function Zl1=Zliq(beta,epsilon,sigma,q)
Zl1=beta+(beta+epsilon.*beta).*(beta+sigma.*beta).*(1./q./beta);
err=1;
while err>eps
    Zl2=beta+(Zl1+epsilon.*beta).*(Zl1+sigma.*beta).*((1+beta-Zl1)./(q.*beta));
    err=(Zl2-Zl1)./Zl2;
    Zl1=Zl2;
end
end
