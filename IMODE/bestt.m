function i = bestt(n,D)
i = zeros(1,n);k = n;
if 2*D > n
D = 1;
n = max(round(0.1*n),2);
end
for j = 1:k
    i(j) = min(randperm(n,D));
end
end
