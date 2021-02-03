function psat=antoine(t)
load a.dat
load b.dat
load c.dat
psat=exp(a-b./(c+t));
end