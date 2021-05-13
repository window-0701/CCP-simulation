function [E,V] = get_V(ni,ne,t)
global n T h beta 
%V = zeros(n+1,1);
av = zeros(n+1,n+1);
bv = zeros(n+1,1);
% V(1) = sin(2 * pi * t);
% V(n) = 0;
%-----------get tridiagonal metrix------------%
av(1,1) = 1;        bv(1) = sin(2 * pi * T * t);
av(n+1,n+1) = 1;    bv(n+1) = 0;
av(2,2) = 2;
av(2,3) = -1;
bv(2) = beta * h * h * (ni(2) - ne(2)) + sin(2 * pi * T * t);

for k = 3:n-1
    av(k,k-1) = -1;
    av(k,k) = 2;
    av(k,k+1) = -1;
    bv(k) =  beta * h * h * (ni(k) - ne(k));
end

av(n,n-1) = -1;
av(n,n) = 2;
bv(n) = beta * h * h * (ni(n) - ne(n));

%------------get new V E---------------%
V = (av\bv);
for k = 1:n
    E(k) = -(V(k+1) - V(k)) / h; 
end

end