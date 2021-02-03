function pexp1=pexp(T)
C1=87.829;C2=-6996.4;C3=-9.8802;C4=7.2099e-6;C5=2;
pexp1=exp(C1+(C2/T)+C3*log(T)+C4*T^C5);
pexp1=pexp1*9.869233e-06;
end