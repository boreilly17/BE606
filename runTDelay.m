% Run tdelay

tspan = [0 30];
lags = [3.3];
history = [0.02; 0.02];

sol = dde23(@tdelaySimple, lags, history, tspan);

figure(1)
plot(sol.x, sol.y)

tspan2 = [0 30];
lags2 = [0.9];
history2 = [0; 0.9; 0];

sol2 = dde23(@modelVI, lags2, history2, tspan2);

figure(2)
plot(sol2.x, sol2.y(1,:), sol2.x, sol2.y(2,:))