function [dxdt] = modelVI_nodelay(x, c)


Bx = c(1);
By = c(2);
ax = c(3);
ay = c(4);
axy = c(5);
as = c(6);
Bs = c(7);
noise = c(8);
n = c(9);
ayd = c(10);
ad = c(11);

dxdt = [
        Bx * (x(3)^n)/((x(3)^n) + 1) *noise - (axy * x(2) * x(1));
        By*x(1) * noise - (ay * x(2)) -  (ayd * x(2) * x(4));
        Bs - as*x(2)*x(3);
        -ad*x(4)
        ];


end