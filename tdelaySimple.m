function [dxdt] = tdelaySimple(t, x, tau)

Bx = 2.3;
By = 24;
ax = 0;
ay = 24;
axy = 120;
noise = 1;

dxdt = [
        Bx*noise - (ax * x(1)) - (axy * x(2) * x(1));
        By*(tau(1,1)) * noise - (ay * x(2))
        ];


end