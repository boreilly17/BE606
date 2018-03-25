%% Normal state of system with signal turned on
% Here, model VI from Geva-Zatorsky et al (2006) is replicated using the
% same paramteters and initial conditions from the paper. We show how the
% model reacts with a time delay when the system has a signal input 

clear all; close all;

% Rate constants determined from paper:
Bx = .9; %Rate of production of p53
By = 1; %Rate of production of mdm2
ax = 0; %Rate of degradation of p53
ay = 0.7; %Rate of degradation of mdm2
axy = 1.4; %Rate of ubiquitination of p53 by mdm2 
            % (assumed immediate degradation if ubiquitinated)
as = 2.7; %Rate of signal degradation
Bs = 0.9; %Signal strength caused by DNA damage
noise = 1; %Small noise variations can be introduced
n = 4; %cooperativity of p53 DNA binding

%Create a constants vector to pass into functions
c = [Bx; By; ax; ay; axy; as; Bs; noise; n];

tspan = [0 24]; %looking at 24 hours of response
lags = 0.9; %This is the time delay for p53 effect on mdm2 production
history = [0; 0.9; 0]; %initial conditions for [p53; mdm2; signal]

%solve the time delayed differential equation over time
solHealthySig = dde23(@modelVI, lags, history, tspan, '', c);

%plots p53 and mdm2
figure(1)
plot(solHealthySig.x, solHealthySig.y(1,:), solHealthySig.x, ...
    solHealthySig.y(2,:), 'LineWidth', 2)
legend('p53', 'mdm2');
title('Typical Stress Response');
xlabel('Time (hours)');
ylabel('Concentrations');



%% Stability analysis of healthy system with signal and no time delay
% In this section we look for steady states of the system, linearize the
% system about those steady states, and determine stability of the steady
% state. Note: at this point, we make the time delay = 0.

%find the jaconbian in symbolic terms so we can linearize the system
syms p m s
dxdt = [
        Bx * (s^n)/((s^n) + 1) *noise - (axy * m * p);
        By*p * noise - (ay * m);
        Bs - as*m*s];
   
J=jacobian(dxdt,[p m s]);
 
% This handy fuction will take our function and find its steady state
% values given the initial conditions
y_steadystate = fsolve(@modelVI_nodelay_noDrug, history, '', c)
 
% Take those values and evaluate the jacobian at steady state
J_eval=  subs(J,[p m s],  y_steadystate');
 
% find the eigen values of this matrix
evals = double(eig(J_eval));
 
disp('The eigenvalues of the healthy system with no delay are: ');
disp(evals);


solHealthySig_noDelay = ode45(@modelVI_nodelay_noDrug_ode, tspan, ...
    history, '', c);

%plots p53 and mdm2
figure(2)
plot(solHealthySig_noDelay.x, solHealthySig_noDelay.y(1,:),...
    solHealthySig_noDelay.x, solHealthySig_noDelay.y(2,:),'LineWidth', 2);
legend('p53', 'mdm2');
title('Typical Stress Response with ''No time-delay'' condition');
xlabel('Time (hours)');
ylabel('Concentrations');



%% Healthy response to dynamic signal
% Here, the same healthy model VI is used from the paper, however this time
% the signal changes. We see how the expression levels adjust in a quick and
% repeatable manner.

flips = 24; %how many times will the signal change
fullSol = [];
timepoints = [];
signalHistory = [];

tspan = [0 5]; %how frequently the signal may change (hours)
history = [0; 0.9; 0];

for i = 1:flips
    if (i > 1 && i < 8) || (i > 10 && i < 14) || (i > 14 && i < 22)
        Bs = 0.9; %signal is on
    else
        Bs = 0; %signal is off
    end
    
    c = [Bx; By; ax; ay; axy; as; Bs; noise; n]; %have to put Bs back
    sol = dde23(@modelVI, lags, history, tspan, '', c);
    
    %store all data in growing vectors
    fullSol = [fullSol, sol.y];
    timepoints = [timepoints, sol.x+((i-1)*tspan(2))];
    signalHistory = [signalHistory, (sol.y*0)+Bs];
    
    history = sol.y(:,end); %set new initial points to end of last run
    
end

figure(3)
subplot(2,1,1)
plot(timepoints, fullSol(1,:), timepoints, fullSol(2,:), 'LineWidth', 2);
legend('p53', 'mdm2');
title('Flipping signal response');
subplot(2,1,2)
plot(timepoints, signalHistory, 'b', 'LineWidth', 2);
legend('Signal');
title('Flipping signal conditions');

%% mdm2 inhibitor

% Change degradation rate of mdm2 such that is overpowers p53


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
plot(timepoints, fullSol(1,:), timepoints, fullSol(2,:),timepoints, fullSol(3,:))
legend('p53', 'mdm2', 'Signal');
title('Tumor Conditions with Drug taken every 24 hrs');

%% STABILITY ANALYSIS CODE - for sick with drug treatment

%So we load in all of our info 
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

%find the jaconbian in symbolic terms
syms p m s d
dxdt = [
        Bx * (s^n)/((s^n) + 1) *noise - (axy * m * p);
        By*p * noise - (ay * m) -  (ayd * m * d);
        Bs - as*m*s;
        -ad*d ];
 J=jacobian(dxdt,[p m s d]);
 
 % This handy fuction will takr our function and find its steady state vals
 y_steadystate=fsolve(@modelVI_nodelay, [0; 0.9; 0;0])
 % Take those values and evaluate the jacobian at steady state
 J_eval=subs(J,[p m s d],  y_steadystate');
 % find the eigen values of this matrix
 [vectors,vals]=eig(J_eval);
 approx_vals=vpa(vals);
STABILITY = [ approx_vals(1,1);approx_vals(2,2); approx_vals(3,3);approx_vals(4,4)]

    % note that these values completely agree with our model
    % we should figure out why the drug doesnt have and imaginary part
    



