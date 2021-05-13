%-------Parameters-------%
global n d T_g p_n dt e V_alpha f gamma n_s T_es T epsilon h beta T_eb_bar D_e_bar D_i_bar D_m_bar mu_e_bar mu_i_bar k_s_bar 
n = 1000; 
d = 2.5; %cm
T_g = 300;%K
p_n = 0.15;%Torr
%--------eV--------------cm3/s-----------eV---%
H_ex = 11.56; k_ex0 = 3.712e-8; E_ex = 15.6;
H_i = 15.7; k_i0 = 1.235e-7; E_i = 18.687;
H_si = 4.14; k_si0 = 2.05e-7; E_si = 4.95;
H_sc = -11.56; k_sc0 = 1.818e-9; E_sc = 2.14;
%----cm3/s-----------cm3/s----------cm3/s---------cm6/s---%
k_r = 2e-7; k_mp = 6.2e-10; k_2q = 3e-15; k_3q = 1.1e-31;
%-----------------        -------------------%
n_nD_e = 3.86e22; n_nD_i = 2.07e18; n_nD_m = 2.42e18;
n_nmu_e = 9.66e21; n_nmu_i = 4.65e19; 
%-----------------        -------------------%
e = 1.6022e-19; epsilon = 8.8542e-14;
n_epsilon = 1e7; T_ei = 1;
gamma = 0.01; k_s = 1.19e7; T_ab = 0.5; V_alpha = 100; f = 13.56; T = 1/f;
n_s = 1e10; T_es = 1; 

dt = 1e-4;%time step
%-----------------        -------------------%
n_n = 4.83e15;
D_e = 7.99171843e6;
D_i = 4.28571429e2;
D_m = 5.010352e2;
mu_e = 2e6;
mu_i = 9.6273292e3;
%-----------------nondimensionalize-------------------%
D_e_bar = 9.42975626e3;
D_i_bar = 5.05689074e-6;
D_m_bar = 5.91191976e-6;
mu_e_bar = 2.35988201;
mu_i_bar = 1.13596805e-2;
%-----------------        -------------------%
alpha = 100;
H_ex_bar = 11.56;
H_i_bar = 15.7;
H_si_bar = 4.14;
H_sc_bar = -11.56;
beta = 1.13096045e3;
n_n_bar = 4.83e5;

E_ex_bar = 15.06;
E_i_bar = 18.687;
E_si_bar = 4.95;
E_sc_bar = 2.14;

k_ex0_bar = 2.73746313e-5;
k_i0_bar = 9.10766962e-5;
k_si0_bar = 1.51179941e-4;
k_sc0_bar = 1.34070796e-6;
k_r_bar = 1.47492625e-4;
k_mp_bar = 4.57227139e-7;
k_2q_bar = 2.21238938e-12;
k_3q_bar = 8.1120944e-19;

n_epsilon_bar = 1e-3;
T_ei_bar = 1;
k_s_bar = 0.351032448;
T_eb_bar = 0.5;
h = 0.001;

%-----------------main-------------------%
%sort: get_ni --> get_ne --> get_nm --> get_p --> get_V

%-----------initial condition------------%
x = 0:h:1;
ne = 16 * n_epsilon_bar * x .* x .* (1 - x).^2;
ni = 16 * n_epsilon_bar * x .* x .* (1 - x).^2;
nm = 16 * n_epsilon_bar * x .* x .* (1 - x).^2;
Te = zeros(n+1,1);
Te = Te + T_ei_bar;
V = zeros(n+1,1);
pe = ne .* Te;

%----------------time step----------------%
for t = 0:dt:1e-1     
    nm(1) = 0;
    nm(n+1) = 0;
    Te(1) = T_eb_bar;
    Te(n+1) = T_eb_bar;
    V(1) = sin(2*pi*t);
    V(n+1) = 0;
    
%-----------------get_V E------------------%
    [E,V] = get_V(ni,ne,t);
    k_ex_bar = n_s ./ f .* k_ex0 .* exp(- E_ex_bar ./ Te);
    k_i_bar = n_s ./ f .* k_i0 .* exp(- E_i_bar ./ Te);
    k_si_bar = n_s ./ f .* k_si0 .* exp(- E_si_bar ./ Te);
    k_sc_bar = n_s ./ f .* k_sc0 .* exp(- E_sc_bar ./Te);
        
    S = k_i_bar .* n_n_bar .* ne + k_si_bar .* nm .* ne + k_mp_bar .* nm .* nm;
    Sm = k_ex_bar .* n_n_bar .* ne - k_si_bar .* nm .* ne - k_sc_bar .* nm .* ne - k_r_bar .* nm .* ne - 2 * k_mp_bar * nm .* nm - k_2q_bar .* n_n_bar .* nm - k_3q_bar .* nm .* n_n_bar .* n_n_bar;
    
    pl = H_ex_bar .* k_ex_bar .* n_n_bar .* ne + H_i_bar .* k_i_bar .* n_n_bar .* ne + H_si_bar .* k_si_bar .* nm .* ne + H_sc_bar .* k_sc_bar .* nm .* ne; 
    
    pa = zeros(n+1,1);
    pa(1) = alpha * E(1) * (gamma * mu_i_bar * ni(1) * E(1) + k_s_bar * ne(1));
    for k = 2:n
       pa(k) = alpha * E(k) * (D_e_bar * (ne(k+1) - ne(k)) / h + mu_e_bar * ne(k) * E(k));
    end
    pa(n+1) = alpha * E(n) * (gamma * mu_i_bar * ni(n+1) * E(n) - k_s_bar * ne(n+1));
    
%-----------------get_ni-------------------%
    ni = get_ni(ni,E,S);
    
%-----------------get_ne-------------------%
    ne = get_ne(ni,ne,S,E);

%-----------------get_nm-------------------%
    nm = get_nm(nm,Sm);
    
%-----------------get_p-------------------%
    [pe,Te] = get_pe(ne,E,pl,pe,pa,Te);
    
end

x = x .* d;
ne = ne .* n_s;
%-----------------plot-------------------%
figure(1);
subplot(2,2,1);
plot(x,ne);
xlabel('x/cm');
ylabel('Electric density/cm-3');
title('ne');
subplot(2,2,2);
plot(x,Te);
title('Te');
xlabel('x/cm');
ylabel('Electric temperature/eV');
subplot(2,2,3);
plot(x,V);
title('V');
xlabel('x/cm');
ylabel('Electric potential/V');
subplot(2,2,4);
plot(x(1:n),E);
title('E');
xlabel('x/cm');
ylabel('Electric field/V cm-1');




