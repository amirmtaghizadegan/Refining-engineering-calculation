function psat=antoine(t)
a=[13.8622 2910.26 216.432];
t=t-273;
psat=exp(a(1)-a(2)/(a(3)+t));
psat=psat/100;
end