%% Clean-up
clear; clear all;

%% Nifty bits
deg2rad = pi/180;

%% System parameters
%Some other constants
V_a = 580/3.6;
aileron_sat = 30*deg2rad; %[deg]
g = 9.81;

%Inner loop transfer function parameters (B&C figure 6.9)
a_2 = -0.65;
a_1 = 2.87;

u_max = 30; %[deg]
e_max = 15; %[deg]

omega_phi = sqrt(abs(a_2)*(u_max/e_max));
seta_phi = 0.707;

k_p_phi = -(u_max/e_max);
k_d_phi = (2*seta_phi*omega_phi - a_1)/a_2;
k_i_phi = -0.1; %Root-locus-curves

omega_chi = omega_phi/20; %Design param
seta_chi = 1; %Design param

k_p_chi = 2*seta_chi*omega_chi*V_a/g;
k_i_chi = (omega_chi^2)*V_a/g;

%% Simulation parameters
sim_time    = 1000;
h           = 0.01;
N           = sim_time/h;
time_vec    = (0:h:sim_time-h);

%% Allocate state vectors
%Specify some initial-values
chi_0       = 15*deg2rad;
p_0         = 0;
phi_0       = 0;

chi         = zeros(1, N); chi(1) = chi_0;
p           = zeros(1, N); p(1) = p_0;
phi         = zeros(1, N); phi(1) = phi_0;
e_chi_int   = zeros(1,N);
e_phi_int   = zeros(1,N);

delta_a     = zeros(1, N);
chi_ref     = zeros(1, N);
chi_ref(N/4:N/2) = 10*deg2rad;
chi_ref(N/2:5*N/8) = 15*deg2rad;
chi_ref(5*N/8:6*N/8) = 5*deg2rad;
d           = 1.5*deg2rad;

%% Simulation loop
for i = 1:N
    t = (i-1)*h;
    e_chi = chi_ref(i) - chi(i);
    phi_c = k_i_chi*e_chi_int(i) + k_p_chi*e_chi;
    e_phi = phi_c - phi(i);
    
    delta_a(i) = k_i_phi*e_phi_int(i) + k_p_phi*e_phi - k_d_phi*p(i);
    if delta_a(i) >= aileron_sat
            delta_a(i) = aileron_sat;
    elseif delta_a(i) <= -aileron_sat
            delta_a(i) = -aileron_sat;
    end
    
    if i < N
        p(i+1)      = p(i) + (delta_a(i)*a_2 - a_1*p(i))*h;
        phi(i+1)    = phi(i) + p(i)*h;
        chi(i+1)    = chi(i) + (g/V_a)*(phi(i) + d)*h;

        e_chi_int(i+1) = e_chi_int(i) + e_chi*h;
        e_phi_int(i+1) = e_phi_int(i) + e_phi*h;
    end
end

%% Plot
figure (1);
hold on;
plot(time_vec, chi_ref, 'b--');
plot(time_vec, chi, 'r');

hold off;
grid on;
title('Course');
xlabel('Time [s]'); 
ylabel('Course [rad]');
legend('Ref', 'Course');

figure(2)
hold on;
plot(time_vec, delta_a, 'r');
hold off;
grid on;
title('Input');
xlabel('Time [s]'); 
ylabel('Aileron [rad]');