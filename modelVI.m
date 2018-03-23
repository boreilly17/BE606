function [dxdt] = modelVI(t, x, tau, c)

Bx = c(1);
By = c(2);
ax = c(3);
ay = c(4);
axy = c(5);
as = c(6);
Bs = c(7);
noise = c(8);
n = c(9);


dxdt = [
        Bx * (x(3)^n)/((x(3)^n) + 1) *noise - (axy * x(2) * x(1));
        By*(tau(1,1)) * noise - (ay * x(2));
        Bs - as*x(2)*x(3)
        ];


end