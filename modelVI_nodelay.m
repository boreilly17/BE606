function [dxdt] = modelVI_nodelay(x, c)


Bx = .9;
By = 1;
ax = 0;
ay = 0.2;
axy = 1.4;
as = 2.7;
Bs = 10;
noise = 1;
n = 4;
Bx = 1 * Bx;
ayd = 10;
ad = 0.5;


dxdt = [
        Bx * (x(3)^n)/((x(3)^n) + 1) *noise - (axy * x(2) * x(1));
        By*x(1) * noise - (ay * x(2)) -  (ayd * x(2) * x(4));
        Bs - as*x(2)*x(3);
        -ad*x(4)
        ];


end