function [pe,Te] = get_pe(ne,E,pl,pe,pa,Te)
global n dt h T_eb_bar D_e_bar mu_e_bar 
%-----------get tridiagonal metrix------------%
ap = zeros(n+1,n+1);
bp = zeros(n+1,1);
ap(1,1) = 1;        ap(n+1,n+1) = 1;
bp(1) = ne(1) * T_eb_bar;
bp(n+1) = ne(n+1) * T_eb_bar;

ap(2,2) = 1 + 5 * dt / (6*h*h) * (4 * D_e_bar - mu_e_bar * h * (E(2) - E(1)));
ap(2,3) = - 5 * dt / (6*h*h) * (2 * D_e_bar + mu_e_bar * h * E(2));
bp(2) = pe(2) + 2 * dt / 3 * (pa(2) - pl(2)) - dt * D_e_bar / (3*h*h) * ((ne(2) + ne(3)) * (Te(3) - Te(2)) - ((ne(1) + ne(2)) * (Te(2) - Te(1)))) + (5 * dt / (6*h*h) * (2 * D_e_bar - mu_e_bar * h * E(1))) * (ne(1) * T_eb_bar);

for k = 3:n-1
    ap(k,k-1) = - 5 * dt / (6*h*h) * (2 * D_e_bar - mu_e_bar * h * E(k-1));
    ap(k,k) = 1 + 5 * dt / (6*h*h) * (4 * D_e_bar - mu_e_bar * h * (E(k) - E(k-1)));
    ap(k,k+1) = - 5 * dt / (6*h*h) * (2 * D_e_bar + mu_e_bar * h * E(k));
    bp(k) = pe(k) + 2 * dt / 3 * (pa(k) - pl(k)) - dt * D_e_bar / (3*h*h) * ((ne(k) + ne(k+1)) * (Te(k+1) - Te(k)) - ((ne(k-1) + ne(k)) * (Te(k) - Te(k-1))));
end
ap(n,n-1) = - 5 * dt / (6*h*h) * (2 * D_e_bar - mu_e_bar * h * E(n-1));
ap(n,n) = 1 + 5 * dt / (6*h*h) * (4 * D_e_bar - mu_e_bar * h * (E(n) - E(n-1)));
bp(n) = pe(n) + 2 * dt / 3 * (pa(n) - pl(n)) - dt * D_e_bar / (3*h*h) * ((ne(n) + ne(n+1)) * (Te(n+1) - Te(n)) - ((ne(n-1) + ne(n)) * (Te(n) - Te(n-1)))) - (- 5 * dt / (6*h*h) * (2 * D_e_bar + mu_e_bar * h * E(n))) * ne(n+1) * T_eb_bar;
%----------get new pe-------------%
pe = ap\bp;
Te = pe./ne;
if Te < 0
    Te = 0;
end

end