%% y steady state 



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

syms p m s d
 
dxdt = [
        Bx * (s^n)/((s^n) + 1) *noise - (axy * m * p);
        By*p * noise - (ay * m) -  (ayd * m * d);
        Bs - as*m*s;
        -ad*d
        ];
 
 
 J=jacobian(eqn1,[p m s d])
 
 yss=fsolve(@function1, [0; 0.9; 0;0])
 
    