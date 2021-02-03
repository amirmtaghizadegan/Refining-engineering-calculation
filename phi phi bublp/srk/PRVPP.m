function psat=PRVPP(T,Tc,Pc,w)
sigma=1;epsilon=0;ohm=0.08664;psi=0.42748;
Tr=T./Tc;
alpha=(1+(0.480+1.574*w-0.176*w.^2).*(1-Tr.^0.5)).^2;
R=82.05;
a=psi.*alpha.*(R.^2.*Tc.^2)./Pc;
b=ohm*R.*Tc./Pc;
psat=antoine(T);
psat=[0 psat];
beta=b.*psat/R/T;
q=a./b./R./T;
Zl=Zliq(beta,epsilon,sigma,q);
Zv=Zvap(beta,q,epsilon,sigma);
Il=1/(sigma-epsilon)*log((Zl+sigma.*beta)./(Zl+epsilon.*beta));
Iv=1/(sigma-epsilon)*log((Zv+sigma.*beta)./(Zv+epsilon.*beta));
lnphil=Zl-1-log(Zl-beta)-q.*Il;
lnphiv=Zv-1-log(Zv-beta)-q.*Iv;
err=1;
while err>eps
    psat=psat.*(exp(lnphil)./exp(lnphiv));
    beta=b.*psat/R/T;
    Zl=Zliq(beta,epsilon,sigma,q);
    Zv=Zvap(beta,q,epsilon,sigma);
    Il=1/(sigma-epsilon)*log((Zl+sigma.*beta)./(Zl+epsilon.*beta));
    Iv=1/(sigma-epsilon)*log((Zv+sigma.*beta)./(Zv+epsilon.*beta));
    lnphil=Zl-1-log(Zl-beta)-q.*Il;
    lnphiv=Zv-1-log(Zv-beta)-q.*Iv;
    err=lnphil-lnphiv;
end
end