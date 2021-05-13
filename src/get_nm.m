function nm = get_nm(nm,Sm)
global n dt h D_m_bar
% nm(1) = 0;
% nm(n) = 0;
%-----------get tridiagonal metrix------------%
am = zeros(n+1,n+1);
bm = zeros(n+1,1);
am(1,1) = 1;	am(n+1,n+1) = 1;
bm(1) = 0;  	bm(n+1) = 0;
am(2,2) = 1 + 2 * D_m_bar * dt / (h * h);
am(2,3) = - D_m_bar * dt / (h*h);
bm(2) = nm(2) + dt * Sm(2);

for k = 3:n-1
    am(k,k-1) = - D_m_bar *dt / (h*h);
    am(k,k) = 1 + 2 * D_m_bar * dt / (h * h);
    am(k,k+1) = - D_m_bar * dt / (h*h);
    bm(k) = nm(k) + dt * Sm(k);
end

am(n,n-1) = - D_m_bar * dt / (h*h);
am(n,n) = 1 + 2 * D_m_bar * dt / (h * h);
bm(n) = nm(n) + dt * Sm(n);

%-------------get new nm---------------%
nm = (am\bm);

end