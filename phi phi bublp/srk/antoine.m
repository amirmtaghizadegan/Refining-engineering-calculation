function psat=antoine(t)
a=13.6608;
b=2154.70;
c=238.789;
t=t-273;
psat=exp((a-b./(c+t)));
psat=psat/100;
end