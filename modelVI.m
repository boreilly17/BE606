function [dxdt] = modelVI(t, x, tau)

Bx = .9;
By = 1;
ax = 0;
ay = .7;
axy = 1.4;
as = 2.7;
Bs = 0.9;
noise = 1;
n = 4;

dxdt = [
        Bx * (x(3)^n)/((x(3)^n) + 1) *noise - (axy * x(2) * x(1));
        By*(tau(1,1)) * noise - (ay * x(2));
        Bs - as*x(2)*x(3)
        ];


end