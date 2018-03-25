%% STABILITY ANALYSIS CODE - for sick with drug treatment

% this is for our code when there is 


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
