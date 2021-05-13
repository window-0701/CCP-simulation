function ni = get_ni(ni,E,S)
global n dt h D_i_bar mu_i_bar
%-----------get tridiagonal metrix------------%
ai = zeros(n+1,n+1);
bi = zeros(n+1,1);
ai(1,1) = 1 + dt / (h*h) * (2 * D_i_bar - h * mu_i_bar * E(1));
ai(1,2) = - dt / (h*h) * (2 * D_i_bar - h * mu_i_bar * E(1));
bi(1) = dt * S(1) + ni(1);
for k = 2:n
    ai(k,k-1) = - dt / (2*h*h) * (2 * D_i_bar + h * mu_i_bar * E(k-1));
    ai(k,k) = 1  + dt / (2*h*h) * (4 * D_i_bar + h * mu_i_bar * (E(k) - E(k-1)));
    ai(k,k+1) = - dt / (2*h*h) * (2 * D_i_bar - h * mu_i_bar * E(k));
    bi(k) = ni(k) + dt * S(k);
end
ai(n+1,n) = - dt / (h*h) * (2 * D_i_bar + h * mu_i_bar * E(n));
ai(n+1,n+1) = 1 + dt / (h*h) * (2 * D_i_bar + h * mu_i_bar * E(n));
bi(n+1) = ni(n+1) + dt * S(n+1);
%------------get new ni----------------------%
ni = (ai\bi);

end