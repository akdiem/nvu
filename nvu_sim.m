clear
clear global

global Jrho_IN x_rel s m cm um  g A mM uM nM C V mV F pF S pS dyn mmHg

% Units
s = 1;
cm = 1e-2; um = 1e-4*cm;
mM = 1e-3; uM = 1e-3 * mM; nM = 1e-6 * mM;
C = 1;
A = C/s;
V = 1; mV = 1e-3 * V;
F = C/V; pF = 1e-12 * F;
S = A/V; pS = 1e-12 * S;
dyn = 1; mmHg = 1333.22 * dyn/cm^2;

% Variables
potassium_s = 2.92655044308714e-8;
%potassium_s = 0;
ip3 = 5.37611796987610e-10;
calcium_a = 1.47220569018281e-07;
h = 0.404507631346124;
ss = 0.0161921297424289;
eet = 4.78801348065449e-07;
nbk = 6.24930194169376e-5;
Vk = -0.0814061063457068;
potassium_p = 0.00353809145071707;
calcium_p = 4.60269585230500e-06;
k = 8.01125818473096e-09;
Vm = 8.33004194103223e-05;
n = 0.283859572906570;
%x = 0.000124217156123987;
x = 2*pi*20*um;
calcium_smc = 3.41385670857693e-07;
omega = 0.536911672725179;
yy = 0.000115089683436595;
y0 = [potassium_s ip3 calcium_a h ss eet nbk Vk potassium_p calcium_p...
    k Vm n x calcium_smc omega yy];

x_rel = x;

% t1=-20;t2=0;
% sizeJrho=160000;
% tspan = [t1,t2]; 
% Max_neural_Kplus = 1.3/100*mM/s;
% Max_neural_glut = 0.75;
% Jrho_IN=zeros(sizeJrho,3);
% Jrho_IN(:,1) = linspace(t1,t2,sizeJrho);
% it1=10000; it2=4000; it3=60000; it4=1000; it5=sizeJrho-(it1+it2+it3+it4);
% Jrho_IN(:,2) = Max_neural_Kplus*[zeros(it1,1); linspace(0,1,it2)'; ones(it3,1); linspace(1,0,it4)';zeros(it5,1)];
% Jrho_IN(:,3) =  Max_neural_glut*[zeros(it1,1); linspace(0,1,it2)'; ones(it3,1); linspace(1,0,it4)';zeros(it5,1)];
% options = odeset('reltol',1e-8,'abstol',1e-8);
% [t, y] = ode23tb(@nvu, tspan, y0, options);
% y0 = y(end,:);


t1=0;t2=40;
sizeJrho=160000;
tspan = [t1,t2]; 
Max_neural_Kplus = 0.48*uM/s;
Max_neural_glut = 0.5;
% Max_neural_Kplus = 0.1*uM/s;
% Max_neural_glut = 0;
Jrho_IN=zeros(sizeJrho,3);
Jrho_IN(:,1) = linspace(t1,t2,sizeJrho);
it1=10000; it2=4000; it3=60000; it4=1000; it5=sizeJrho-(it1+it2+it3+it4);
Jrho_IN(:,2) = Max_neural_Kplus*[zeros(it1,1); linspace(0,1,it2)';...
    ones(it3,1); linspace(1,0,it4)';zeros(it5,1)];
Jrho_IN(:,3) =  Max_neural_glut*[zeros(it1,1); linspace(0,1,it2)';...
    ones(it3,1); linspace(1,0,it4)';zeros(it5,1)];

% Stiff solver
options = odeset('reltol',1e-8,'abstol',1e-8);
[t, y] = ode23tb(@nvu, tspan, y0, options);

figure(1);
subplot(4,2,1);
plot(t,y(:,1)./uM)
ylabel("K+ syn [uM]");
subplot(4,2,3);
plot(t,y(:,2)./uM)
ylabel("IP3 [uM]");
subplot(4,2,5);
plot(t,y(:,3)./uM)
ylabel("Ca2+ ast [uM]");
subplot(4,2,7);
plot(t,y(:,6)./uM)
ylabel("EET [uM]");
subplot(4,2,2);
plot(t,y(:,8)./mV)
ylabel("Vk [mV]");
subplot(4,2,4);
plot(t,y(:,9)./mM)
ylabel("K+ ast [mM]");
subplot(4,2,6);
plot(t,y(:,15)./uM)
ylabel("Ca2+ smc [uM]");
subplot(4,2,8);
plot(t,y(:,14)./(2*pi)./um)
ylabel("radius [um]");


function dydt = nvu(t, y)
% This function calculates ion currents in the NVU

global Jrho_IN x_rel s m cm um  g A mM uM nM C V mV F pF S pS dyn mmHg

% Variables
%potassium_s = K_release(t)*uM;
potassium_s = y(1);
ip3 = y(2);
calcium_a = y(3);
h = y(4);
ss = y(5);
eet = y(6);
nbk = y(7);
Vk = y(8);
potassium_p = y(9);
calcium_p = y(10);
k = y(11);
Vm = y(12);
n = y(13);
x = y(14);
calcium_smc = y(15);
omega = y(16);
yy = y(17);

% Parameter
JSigKkNa = 0.6 * mM/s;
KKoa = 1.5 * mM;
rh = 4.8 * uM;
kdeg = 1.25 * 1/s;
KG = 8.82;
delta = 0.001235;
beta = 0.0244;
Jmax = 2880 * uM/s;
Ki = 0.03 * uM;
Kact = 0.17 * uM;
Vmax = 20 * uM/s;
Kp = 0.24 * uM;
%Kp = 0.192 * uM;
Pl = 5.2 * uM/s;
calcium_er = 400 * uM; % value from Bennet et al.
Castr = 40 * pF;
gamma = 1970 * mV/uM;
%gamma = 834.3/2 * mV/uM;
gtrpv = 200 * pS;
% gtrpv = 50 * pS;
vtrpv = 6 * mV;
kon = 2 * 1/(uM*s);
Kinh = 0.1 * uM;
tautrpv = 0.9 * 1/s;
eps12 = 0.16;
%eps12 = 0.1;
kappa = 0.04;
%kappa = 0.1;
gammacai = 0.2 * uM;
%gammacai = 0.01 * uM;
gammacae = 0.2 * mM;
v1trpv = 120 * mV;
v2trpv = 13 * mV;
Veet = 72 * 1/s;
calcium_a_min = 0.1 * uM;
keet = 7.1 * 1/s;
psibk = 2.664 * 1/s;
v4bk = 14.5 * mV;
v5bk = 8 * mV;
v6bk = -15 * mV;
eetshift = 2 * mV/uM;
Ca3bk = 400 * nM;
Ca4bk = 150 * nM;
gbk = 225.6 * pS;
vbk = -95 * mV;
%vbk = -70 * mV;
gleak = 78.54 * pS;
vleak = -70 * mV;
%vleak = -60 * mV;
Csmc = 19.635 * pF;
gkir0 = 145 * pS;
vkir1 = 57 * mV;
vkir2 = 130 * mV;
VRpa = 3.2e-5;
VRps = 0.1;
gca = 157 * pS;
vca = 80 * mV;
v2 = 25 * mV;
dp = 60 * mmHg;
Rdecay = 1 * 1/s;
potassium_p_min = 3 * mM;
Cadecay = 0.5 * 1/s;
calcium_p_min = 2 * uM;
alphakir = 1020 * s;
av1 = 18 * mV;
av2 = 10.8 * mV;
betakir = 26.9 * s;
bv1 = 18 * mV;
bv2 = 0.06 * mV;
gl = 62.832 * pS;
vl = -70 * mV;
gk = 251.33 * pS;
vk = -80 * mV;
phin = 2.664;
Ca3 = 400 * nM;
Ca4 = 150 * nM;
v4 = 14.5 * mV;
v5 = 8 * mV;
v6 = -15 * mV;
tau = 0.2 * dyn/cm;
Kd = 1000 * nM;
Bt = 10000 * nM;
alpha = 4.3987e15 * nM/C;
kca = 1.3568e11 * nM/C;
we = 0.9;
x0 = 188.5 * um;
x1 = 1.2;
x2 = 0.13;
x3 = 2.2443;
x4 = 0.71182;
x5 = 0.8;
x6 = 0.01;
x7 = 0.32134;
x8 = 0.88977;
x9 = 0.0090463;
sigma0h = 3e6 * dyn/cm^2;
u1 = 41.76;
u2 = 0.047396;
u3 = 0.0584;
wm = 0.7;
kpsi = 3.3;
Cam = 500 * nM;
psim = 0.3;
q = 3;
y0 = 0.928;
y1 = 0.639;
y2 = 0.35;
y3 = 0.78847;
y4 = 0.8;
sigmay0h = 2.6e6 * dyn/cm^2;
vref = 0.24;
Caref = 510 * nM;
ad = 0.28125;
bd = 5;
cd = 0.03;
dd = 1.3;
Sx = 40000 * um^2;
Ax = (x_rel/(2*pi))^2 * pi;

% Synaptic space
JSigK = JSigKkNa * potassium_s/(potassium_s + KKoa);
JKss = interp1(Jrho_IN(:,1),Jrho_IN(:,2),t);
potassium_s_dt = JKss - JSigK;
%potassium_s_dt = 0;


% Astrocytic space
rhos = interp1(Jrho_IN(:,1),Jrho_IN(:,3),t);
G = (rhos+delta)/(KG + rhos + delta);
ip3_dt = rh*G - kdeg*ip3;
Jip3 = Jmax * ((ip3/(ip3+Ki)) * (calcium_a/(calcium_a+Kact)) * h)^3 *...
    (1 - calcium_a/calcium_er);
Jpump = Vmax * calcium_a^2 / (calcium_a^2 + Kp^2);
Jleak = Pl * (1 - calcium_a/calcium_er);
Itrpv = gtrpv * ss * (Vk - vtrpv);
Jtrpv = -Itrpv/(Castr*gamma);
calcium_a_dt = beta * (Jip3 - Jpump + Jleak + Jtrpv);
h_dt = kon * (Kinh - (calcium_a + Kinh) * h);
tauca = tautrpv * (calcium_p/mM);
eps = (x - x_rel)/x_rel;
Hca = calcium_a/gammacai + calcium_p/gammacae;
sinf = (1/(1 + exp(-(eps-eps12)/kappa))) * (1/(1+Hca) *...
    (Hca + tanh((Vk - v1trpv)/v2trpv)));
ss_dt = 1/tauca * (sinf - ss);
eet_dt = Veet * (calcium_a - calcium_a_min) - keet*eet;
v3bk = -(v5bk/2) * tanh((calcium_a-Ca3bk)/Ca4bk) + v6bk;
phibk = psibk * cosh((Vk-v3bk)/(2*v4bk));
ninf = 0.5 * (1 + tanh((Vk + eetshift*eet - v3bk)/v4bk));
nbk_dt = phibk * (ninf - nbk);
Ibk = gbk * nbk * (Vk - vbk);
Ileak = gleak * (Vk - vleak);
Isigk = -JSigK * Castr * gamma;
Vk_dt = 1/Castr * (-Isigk - Ibk - Ileak - Itrpv);


% Perivascular space
Jbk = Ibk/(Castr*gamma);
gkir = gkir0 * sqrt(potassium_p/mM);
vkir = vkir1 * log10(potassium_p/mM) - vkir2;
Ikir = gkir * k * (Vm - vkir);
Jkir = Ikir/(Csmc*gamma);
potassium_p_dt = Jbk/VRpa + Jkir/VRps - Rdecay * (potassium_p -...
    potassium_p_min);
v1 = (-17.4-12*(dp/mmHg)/200)*mV;
minf = 0.5 * (1 + tanh((Vm-v1)/v2));
Ica = gca * minf * (Vm - vca);
Jca = -Ica/(Csmc*gamma);
calcium_p_dt = -Jtrpv - Jca - Cadecay * (calcium_p - calcium_p_min);


% Ion currents
alphak = alphakir / (1 + exp((Vm - vkir + av1)/av2));
betak = betakir * exp(bv2/mV * (Vm - vkir + bv1)/mV);
tauk = 1/(alphak+betak);
kinf = alphak/(alphak+betak);
k_dt = 1/tauk * (kinf - k);
Il = gl * (Vm - vl);
Ik = gk * n * (Vm - vk);
Vm_dt = 1/Csmc * (-Il - Ik - Ica - Ikir);
v3 = -(v5/2) * tanh((calcium_smc-Ca3)/Ca4) + v6;
lamn = phin * cosh(0.5*(Vm-v3)/v4);
ninf = 0.5 * (1 + tanh((Vm-v3)/v4));
n_dt = lamn * (ninf - n);

% Vessel SMC calcium
rho_smc = (Kd+calcium_smc)^2/((Kd+calcium_smc)^2 + Kd*Bt);
calcium_smc_dt = -rho_smc * (Ica*alpha + kca*calcium_smc);


% Vessel mechanics
fdp = 0.5 * dp * (x/pi - Ax/x) * cm;
xd = x/x0;
sigmax = x3*(1 + tanh((xd-x1)/x2)) + x4*(xd-x5) - x8*(x6/(xd-x7))^2 - x9;
fx = we*Sx*sigmax*sigma0h;
u = x-yy;
ud = u/x0;
sigmau = u2 * exp(u1*ud) - u3;
fu = wm*Sx*sigmau*sigma0h;
x_dt = 1/tau * (fdp - fx - fu);
%x_dt = 0;
psi = calcium_smc^q/(Cam^q+calcium_smc^q);
omega_dt = kpsi * (psi/(psim+psi) - omega);
yd = yy/x0;
psiref = Caref^q/(Cam^q+Caref^q);
omega_ref = psiref/(psim + psiref);
sigmay0 = sigmay0h * omega/omega_ref;
sy = (y1/(yd+y2))^y4;
sigmay = sigmay0/sigma0h * (exp(-(yd-y0)^2/(2*sy^2)) - y3)/(1-y3);
cond = sigmau/sigmay;
if cond < 1
    ycond = -vref * psi/psiref * ad * (1-cond)/(ad+cond);
else
    ycond = cd*(exp(bd * (cond-dd)) - exp(bd * (1-dd)));
end
yy_dt = x0*ycond;

dydt = [potassium_s_dt; ip3_dt; calcium_a_dt; h_dt; ss_dt; eet_dt;...
    nbk_dt; Vk_dt; potassium_p_dt; calcium_p_dt; k_dt; Vm_dt; n_dt;...
    x_dt; calcium_smc_dt; omega_dt; yy_dt];
end


function c = K_release(t)
if t <= 20
    c = 1.4 * (1+tanh(t/2-2))/2;
else
    c = 1.4 * exp(-(t-20)/2);
end
end


function c = Glu_release(t)
if t <= 20
    c = 0.2 * (1+tanh(t/2-2))/2;
else
    c = 0.2 * exp(-(t-20)/2);
end
end