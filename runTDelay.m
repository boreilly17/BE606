% Run tdelay


%% Typical Stress Response
%Constants

Bx = .9;
By = 1;
ax = 0;
ay = 0.7;
axy = 1.4;
as = 2.7;
Bs = 0.9;
noise = 1;
n = 4;

c = [Bx; By; ax; ay; axy; as; Bs; noise; n];

tspan = [0 30];
lags = [0.9];
history = [0; 0; 0];

sol = dde23(@modelVI, lags, history, tspan, '', c);

figure(1)
plot(sol.x, sol.y(1,:), sol.x, sol.y(2,:))
legend('p53', 'mdm2');
title('Typical Stress Response');

%% Tumor onset

% Change degradation rate of mdm2 such that is overpowers p53
clear 

tspan = [0 30];
lags = [0.9];
history = [0; 0.9; 0];

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
c = [Bx; By; ax; ay; axy; as; Bs; noise; n];
sol = dde23(@modelVI, lags, history, tspan, '', c);

figure(2)
plot(sol.x, sol.y(1,:), sol.x, sol.y(2,:))
legend('p53', 'mdm2');
title('Tumor Conditions');



%% mdm2 inhibitor

% Change degradation rate of mdm2 such that is overpowers p53
clear 

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


c = [Bx; By; ax; ay; axy; as; Bs; noise; n; ayd; ad];

flips = 4;
fullSol = [];
timepoints = [];

tspan = [0 24];
lags = [0.9];
history = [0; 0.9; 0; 0];

for i = 1:flips
    if i == 1
        drugDose = 0;
    else
        drugDose = 100;
    end
    
    sol = dde23(@modelVI_drug, lags, history, tspan, '', c);
    
    fullSol = [fullSol, sol.y];
    timepoints = [timepoints, sol.x+((i-1)*tspan(2))];
    
    history = sol.y(:,end);
    history(4) = 100;
    
end


figure(2)
plot(timepoints, fullSol(1,:), timepoints, fullSol(2,:))
legend('p53', 'mdm2');
title('Tumor Conditions with Drug taken every 24 hrs');


%% Healthy changes in signal

flips = 4;
fullSol = [];
timepoints = [];

tspan = [0 60];

for i = 1:flips
    if i == 1
        Bs = 0.1;
    elseif i ==2
        Bs = .5;
    elseif i ==3
        Bs = .7;
    else
        Bs = .1;
    end
    
    c = [Bx; By; ax; ay; axy; as; Bs; noise; n];
    sol = dde23(@modelVI, lags, history, tspan, '', c);
    
    fullSol = [fullSol, sol.y];
    timepoints = [timepoints, sol.x+((i-1)*tspan(2))];
    
    history = sol.y(:,end);
    
end



figure(3)
plot(timepoints, fullSol(1,:), timepoints, fullSol(2,:))
legend('p53', 'mdm2');
title('Flipping signal conditions');

