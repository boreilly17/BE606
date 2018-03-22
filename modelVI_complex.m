function [dxdt] = modelVI_complex(t, x, tau)

theta_damage = 1;

n = 4;

T_s = 1;
B_p53 = 0.9;
B_inhibitor = 1;
B_mdm2 = B_inhibitor;
B_signal = 0.9;
alpha_signal = 2.7;
alpha_inhibitor = 0.7;
alpha_mdm2 = alpha_inhibitor;
alpha_active = 1.4;
alpha_inactive = 100*alpha_active;
w = 500 * alpha_active;

p53_lag = tau(1);

p53_inactive = x(1);
p53_active = x(2);
mdm2 = x(3);
I = x(4);
S = x(5);

%{
species list:
p53_inactive
p53_active
mdm2
I (inhibitor)
S (signal)
%}

dxdt = [
        B_p53 - (alpha_inactive * mdm2 * (p53_active) - ...
            (w * ((S^n)/((S^n) + T_s)) * (p53_active)));
        (w * ((S^n)/((S^n) + T_s)) * p53_inactive) - ...
            (alpha_active * mdm2 * p53_active);
        (B_mdm2 * p53_lag) - (alpha_mdm2 * mdm2);
        (B_inhibitor * p53_lag) - (alpha_inhibitor * I);
        (B_signal * theta_damage) - (alpha_signal * I * S)
        ];


end