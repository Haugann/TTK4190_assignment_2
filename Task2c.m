clear; clc;
%% Variable declarations
a1_phi = 2.87;
a2_phi = -0.65;

d_a_max = 30;
e_phi_max = 15;
zeta_phi = 0.707;
omega_phi = sqrt(2 * abs(a2_phi));

Kp_phi = d_a_max / e_phi_max * sign(a2_phi);
Kd_phi = (2 * zeta_phi * omega_phi- a1_phi) / a2_phi;

%% Transfer function
H = tf(a2_phi, [1 (a1_phi + Kd_phi * a2_phi) (a2_phi * Kp_phi) 0]);

K_i = linspace(0,10,10000); 
[R, K] = rlocus(H, -K_i);

%% Finding the max gain
for i = 1:length(K)
    if (real(R(1,i)) > 0) || (real(R(2,i)) > 0) || (real(R(3,i)) > 0)
        K_i_max = K(i);
        break
    end
end

%% Generating root locus plot
figure()

subplot(1,2,1)
rlocusplot(H, K_i, 'b');
title('Positive values of $K_{\phi_i}$','Interpreter','latex', 'fontsize', 16);
xlabel('Im'); ylabel('Re')

subplot(1,2,2)
rlocusplot(H, -K_i, 'r');
title('Negative values of $K_{\phi_i}$','Interpreter','latex', 'fontsize', 16);
xlabel('Im'); ylabel('Re')

set(findall(gcf,'Type','line'),'LineWidth', 1.5)