function ne = get_ne(S,ni,ne,E)
global n dt gamma h D_e_bar mu_e_bar mu_i_bar k_s_bar
%-----------get tridiagonal metrix------------%
ae = zeros(n+1,n+1);
be = zeros(n+1,1);
ae(1,1) = 1 + dt / (h*h) * (2 * D_e_bar + 2 * h * k_s_bar - h * mu_e_bar * E(1));
ae(1,2) = -dt / (h*h) * (2 * D_e_bar +  mu_e_bar * h * E(1));
be(1) = ne(1) + dt * (S(1) - 2 * gamma * mu_i_bar * ni(1) * E(1) / h);
for k = 2:n
    ae(k,k-1) = - dt / (2*h*h) * (2 * D_e_bar - mu_e_bar * h * E(k-1));
    ae(k,k) = 1 + dt / (2 * h*h) * (4 * D_e_bar - mu_e_bar * h *(E(k) - E(k-1)));
    ae(k,k+1) = - dt / (2*h*h) * (2 * D_e_bar + mu_e_bar * h * E(k));
    be(k) = ne(k) + dt * S(k);
end
ae(n+1,n) = - dt / (h*h) * (2 * D_e_bar - mu_e_bar * h * E(n));
ae(n+1,n+1) = 1 + dt / (h*h) * (2 * D_e_bar + 2 * h * k_s_bar + h * mu_e_bar * E(n));
be(n+1) = ne(n+1) + dt * (S(n+1) + 2 * gamma / h * mu_i_bar * E(n) * ni(n+1));
%------------------get new ne-----------------%
ne = (ae\be);

end