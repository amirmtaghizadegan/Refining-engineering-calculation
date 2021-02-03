function virialdat(n)
B=zeros(n);
for i=1:n
    for j=1:n
        fprintf('insert B%d%d (virial coefs.)',i,j)
        B(i,j)=input('');
    end
end
save Bv.dat B -ascii
end