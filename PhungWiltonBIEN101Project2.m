%% *FIGURE 1*

% Define initial conditions and rate constants 
initialConditions =[5e-12 9e-8 1.7e-7 2e-8 7e-10 0.0000014 0 0 0 0 0 0 0 0 0 0 0 0]; % starting conditions (concentrations) 
k = [2e7 2e7 1e7 2e6 1e7 1e8 1e7 4e8 .0005 .4 .3 1.15 8.2 32 1e5 24 44 .001 70 .02]; % rate constants 
time = [0 250];
opts = odeset('RelTol',1e-12,'AbsTol',1e-12); 
[t,y] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);
thrombin = (y(:,9) + 1.2 * y(:,11)) / 1.4e-6 * 100;

% Figure 1A
thrombinA = thrombin;
f1 = figure;
subplot(2,1,1)
plot(t, thrombinA);
xlabel("Time (s)");
ylabel("% Thrombin Formation");
legend("Basic Model", "Location", "Best");
title("Figure 1A");

% Figure 1B 
thrombinB1 = (thrombinA/100) * 1.4; % Thrombin concentration ( not as a percentage like previous ) 

k(7) = 1e6; % change rate constant ( stated in legend ) 
k(9) = 0.0005; % change rate constant ( stated in legend ) 

[t1b, y1b] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts) ;
thrombinB2 = (y1b(:,9) + 1.2 * y1b(:,11)) * 1e6; % in \muM
 
k(7) = 1e7; % reset to original rate constant value
k(9) = 0.05; % reset to original rate constant value

k(8) = 4e7; % change rate constant ( stated in legend ) 
k(10) = 0.04; % change rate constant ( stated in legend )

[t1c, y1c] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);
thrombinB3 = (y1c(:,9) + 1.2 * y1c(:,11)) * 1e6; % in \muM
 
subplot(2,1,2)
plot(t, thrombinB1);
hold on 
plot(t1b, thrombinB2);
plot(t1c, thrombinB3);
legend("Basic Model", "k7 = 1e6, k9 = 0.0005", "k8 = 4e7, k10 = 0.04" , "Location", "Best") ;
xlabel("Time (s)") ;
ylabel("Thrombin Formation (\muM)");
title("Figure 1B"); 
hold off 
saveas(f1, "Figure1.png");

%% *FIGURE 2*

% Reset parameters by redefining 
initialConditions =[ 5e-12 9e-8 1.7e-7 2e-8 7e-10 1.4e-6 0 0 0 0 0 0 0 0 0 0 0 0 ];
k = [ 2e7 2e7 1e7 2e6 1e7 1e8 1e7 4e8 5e-4 4e-1 3e-1 0.5 8.2 32 1e5 24 44 1e-3 70 .02 ]; 
time = [0 250];
opts = odeset('RelTol',1e-12,'AbsTol',1e-12); 
[t,y] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);
k(20) = 0;
[t2,y2] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);

thrombin = (y(:,9) + 1.2 * y(:,11)) / 1.4e-6 * 100;
thrombin2 = (y2(:,9) + 1.2 * y2(:,11)) / 1.4e-6 * 100; 

factorXa = (y(:,16) + y(:,8) + y(:,10)) / 170e-9 * 100; 
factorXa2 = (y2(:,16) + y2(:,8) + y2(:,10)) / 170e-9 * 100;

factorIXa = (y(:,15) + y(:,14) + y(:,7)) / 90e-9 * 100;
factorIXa2 = (y2(:,15) + y2(:,14) + y2(:,7)) / 90e-9 * 100;

f2 = figure; % Decedided to change code to display all figures on same figure but in separare graphs (subplots)

% Figure 2A - Thrombin formation
subplot(3, 1, 1);
plot(t, thrombin);
hold on;
plot(t2, thrombin2);
legend("Basic Model","Alternate Model", "Location", "Best");
xlabel("Time (s)");
ylabel("% Thrombin formation");
title("Figure 2A");
hold off 

% Figure 2B - Factor Xa formation
subplot(3, 1, 2);
plot(t, factorXa);
hold on;
plot(t2, factorXa2);
legend("Basic Model","Alternate Model", "Location", "Best");
xlabel("Time (s)");
ylabel("% Factor Xa formation");
title("Figure 2B");
hold off

% Figure 2C - % Factor IXa formation
subplot(3, 1, 3);
plot(t, factorIXa);
hold on;
plot(t2, factorIXa2);
legend("Basic Model","Alternate Model", "Location", "Best");
xlabel("Time (s)");
ylabel("% Factor IXa formation");
title("Figure 2C");
hold off 
saveas(f2,"Figure2.png");

%% *FIGURE 3*

% Reset parameters by redefining 
initialConditions =[ 5e-12 9e-8 1.7e-7 2e-8 7e-10 1.4e-6 0 0 0 0 0 0 0 0 0 0 0 0 ];
k = [ 2e7 2e7 1e7 2e6 1e7 1e8 1e7 4e8 5e-4 4e-1 3e-1 0.5 8.2 32 1e5 24 44 1e-3 70 .02 ]; 
time = [0 250];
opts = odeset('RelTol',1e-12,'AbsTol',1e-12); 
[t,y] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);
[t3,y3] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);

thrombin3 = (y3(:,9) + 1.2 * y3(:,11)) / 1.4e-6 * 100; 
factorXa3 = (y3(:,16) + y3(:,8) + y3(:,10)) / 170e-9 * 100;
factorVIIIa3 = (y3(:,15) + y3(:,14) + y3(:,7)) / 0.7e-9 * 100;
factorIXa3 = (y3(:,15) + y3(:,14) + y3(:,7)) / 90e-9 * 100;
factorV3 = (y3(:,4)) / 20e-9 * 100 ; 
factorVIII3 = (y3(:,5)) / 0.7e-9 * 100; 

% Figure 3 
f3 = figure; 
plot(t3, thrombin3)
hold on
plot(t3, factorVIII3)
plot(t3, factorV3)
plot(t3, factorXa3)
plot(t3, factorIXa3)
legend("Thrombin","FactorVIII", "FactorV", "FactorXa", "FactorIXa", "Location", "Best")
xlabel("Time (s)");
ylabel("% Formation / Degredation");
title("Figure 3")
hold off
saveas(f3,"Figure3.png");

%% *FIGURE 4* 

% Reset parameters by redefining 
initialConditions =[ 5e-12 9e-8 1.7e-7 2e-8 7e-10 1.4e-6 0 0 0 0 0 0 0 0 0 0 0 0 ];
k = [ 2e7 2e7 1e7 2e6 1e7 1e8 1e7 4e8 5e-4 4e-1 3e-1 0.5 8.2 32 1e5 24 44 1e-3 70 .02 ]; 
time = [0 250];
opts = odeset('RelTol',1e-12,'AbsTol',1e-12); 
[t,y] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);

% [TFVIIa] = 10 pM (10e-12 M)
initialConditions(1) = 10e-12 ; % Changes value for initial concentration
[t4a, y4a] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);
thrombin4a = (y4a(:,9) + 1.2 * y4a(:,11)) * 1e6;

% [TFVIIal = 50 pM (50e-12 M)
initialConditions(1) = 50e-12 ;
[t4b, y4b] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts); 
thrombin4b = (y4b(:,9) + 1.2 * y4b(:, 11)) * 1e6 ;

% [TFVIIal = 500 pM (500e-12 M)
initialConditions(1) = 500e-12;
[t4c, y4c] = ode23s(@(t,y) coagulation(t, y, k),time, initialConditions, opts); 
thrombin4c = (y4c(:,9) + 1.2 * y4c(:,11)) * 1e6; 

% [TFVIIal = 5 nM (5e-9 M)
initialConditions(1) = 5e-9;
[t4d, y4d] = ode23s(@(t,y) coagulation(t, y, k),time, initialConditions, opts); 
thrombin4d = (y4d(:,9) + 1.2 * y4d(:,11)) * 1e6; 

thrombin = thrombin / 100 * 1.4; % Thrombin concentration in \muM 

% Figure 4 
f4 = figure;
plot(t, thrombin);
hold on 
plot(t4a, thrombin4a); 
plot(t4b, thrombin4b);
plot(t4c, thrombin4c);
plot(t4d, thrombin4d);
legend("[TFVIIa] = 5 pM", "[TFVIIa] = 10 pM", "[TFVIIa] = 50 pM", "[TFVIIa] = 500 pM", "[TFVIIa] = 5 nM", 'Location', 'Best');
xlabel("Time (s)");
ylabel("Thrombin Formation (\muM) ");
title("Figure 4");
hold off
saveas(f4, "Figure4.png");

%% *FIGURE 5* 

% Reset parameters by redefining 
initialConditions =[ 5e-12 9e-8 1.7e-7 2e-8 7e-10 1.4e-6 0 0 0 0 0 0 0 0 0 0 0 0 ];
k = [ 2e7 2e7 1e7 2e6 1e7 1e8 1e7 4e8 5e-4 4e-1 3e-1 0.5 8.2 32 1e5 24 44 1e-3 70 .02 ]; 
time = [0 250];
opts = odeset('RelTol',1e-12,'AbsTol',1e-12); 
[t,y] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);

% [VIII] = 0
initialConditions(1) = 5e-12;
initialConditions(5) = 0;
[t5a1, y5a1] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);
thrombin5a1 = (y5a1(:,9) + 1.2 * y5a1(:,11)) * 1e6; 

% [V] = 0
initialConditions(5) = 0.7e-9;
initialConditions(4) = 0;
[t5a2, y5a2] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);
thrombin5a2 = (y5a2(:,9) + 1.2 * y5a2(:,11)) * 1e6; thrombin = (y(:,9) + 1.2 * y(:,11)) / 1e-6;

initialConditions(4) = 20e-9; % Reset to the original value of [V]

% [TFVIIa] = 1nM
initialConditions(1) = 1e-9;
[t5b1, y5b1] =  ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);
thrombin5b = (y5b1(:,9) + 1.2 * y5b1(:,11))*1e6;

% [VIII] = 0
initialConditions(5) = 0;
[t5b2, y5b2] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);
thrombin5b2 = (y5b2(:,9) + 1.2 * y5b2(:,11))*1e6;

% [V] = 0
initialConditions(5) = 0.7e-9;
initialConditions(4) = 0;
[t5b_3, y5b_3] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);
thrombin5b_3 = (y5b_3(:,9) + 1.2 * y5b_3(:,11))*1e6;

f5 = figure;

% Figure 5A

subplot(2,1,1)
plot(t, thrombin);
hold on
plot(t5a1, thrombin5a1);
plot(t5a2, thrombin5a2);
legend("Basic Model", "[VIII] = 0 M","[V] = 0 M", "Location", "Best");
xlabel("Time (s)");
ylabel("Thrombin Formation (\muM) ");
title("Figure 5A");
hold off 

% Figure 5B
subplot(2,1,2)
plot(t5b1, thrombin5b);
hold on
plot(t5b2, thrombin5b2);
plot(t5b_3, thrombin5b_3);
legend("[TFVIIa] = 1 nM", "[VIII] = 0 M", "[V] = 0 M" , "Location", "Best");
xlabel("Time (s)");
ylabel("Thrombin Formation (uM) ");
title("Figure 5B");
hold off
saveas(f5, "Figure5.png");

%% *FIGURE 6*

% Reset parameters by redefining 
initialConditions =[ 5e-12 9e-8 1.7e-7 2e-8 7e-10 1.4e-6 0 0 0 0 0 0 0 0 0 0 0 0 ];
k = [ 2e7 2e7 1e7 2e6 1e7 1e8 1e7 4e8 5e-4 4e-1 3e-1 0.5 8.2 32 1e5 24 44 1e-3 70 .02 ]; 
time = [0 250];
opts = odeset('RelTol',1e-12,'AbsTol',1e-12); 
[t,y] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);

thrombin = (y(:,9) + 1.2 * y(:,11)) / 1e-6;

% k2 = 0
k(2) = 0; % set to new value, 0
[t6a, y6a] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);
thrombin6a = (y6a(:,9) + 1.2 * y6a(:,11))*1e6;

% k1 = 0
k(2) = 2e7; % reset back to original value
k(1) = 0; % set to new value, 0
[t6b, y6b] = ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);
thrombin6b = (y6b(:,9) + 1.2 * y6b(:,11))*1e6;

% k4 = 0
k(1) = 2e7; % reset back to original value
k(4) = 0; % set to new value, 0
[t6c, y6c] =  ode23s(@(t,y) coagulation(t, y, k), time, initialConditions, opts);
thrombin6c = (y6c(:,9) + 1.2 * y6c(:,11))*1e6;

% Figure 6 - Thrombin Formation 
f6 = figure;
plot(t, thrombin);
hold on
plot(t6a, thrombin6a);
plot(t6b, thrombin6b);
plot(t6c, thrombin6c);
legend("Basic Model", "k2 = 0", "k1 = 0", "k4 = 0" ,"Location", "Best");
xlabel("Time (s)");
ylabel("Thrombin Formation (\muM) ");
title("Figure 6");
saveas(f6, "Figure6.png");
hold off

%% Collaborators : Erick Garia, Nohemi Rodriguez Hernandez 
%% Coagulation 
function ydot = coagulation(t, y, k)
ydot = zeros(18,1);

k1=k(1);
k2=k(2);
k3=k(3);
k4=k(4);
k5=k(5);
k6=k(6);
k7=k(7);
k8=k(8);
k9=k(9);
k10=k(10);
k11=k(11);
k12=k(12);
k13=k(13);
k14=k(14);
k15=k(15);
k16=k(16);
k17=k(17);
k18=k(18);
k19=k(19);
k20=k(20);

TFVIIa = y(1);
IX = y(2);
X = y(3);
V = y(4);
VIII = y(5);
II = y(6);
VIIIaIXa = y(7);
VaXa = y(8);
IIa = y(9);
VaXaII = y(10);
mIIa = y(11);
TFVIIaIX = y(12);
TFVIIaX = y(13);
VIIIaIXaX = y(14);
IXa = y(15);
Xa = y(16);
Va = y(17);
VIIIa = y(18);

ydot(1)=k11 * TFVIIaIX - k6 * TFVIIa * IX + k16 * TFVIIaIX + k12 ...
    * TFVIIaX - k6 * TFVIIa * X + k17 * TFVIIaX;

ydot(2)=k16*TFVIIaIX - k6*TFVIIa*IX - k15*IX*Xa - k15*IX*VaXa;

ydot(3)=k17*TFVIIaX - k6*TFVIIa*X - k6*VIIIaIXa*X + k18*VIIIaIXaX;

ydot(4)=-k1*V*Xa - k2*V*IIa - k2*V*mIIa;

ydot(5)=-k3*VIII*Xa - k4*VIII*IIa - k4*VIII*mIIa;

ydot(6)=k19*VaXaII - k6*VaXa*II;

ydot(7)= k7 * VIIIa * IXa - k9 * VIIIaIXa - k6 * VIIIaIXa * X + k18 ...
    * VIIIaIXaX + k13 * VIIIaIXaX - k20 * VIIIaIXa;

ydot(8)=k8*Xa*Va - k10*VaXa + k19*VaXaII - k6*VaXa*II + k14*VaXaII;

ydot(9)=k5*VaXa*mIIa;

ydot(10)=k6*VaXa*II - k19*VaXaII - k14*VaXaII;

ydot(11)=k14*VaXaII - k5*VaXa*mIIa;

ydot(12)=k6*TFVIIa*IX - k16*TFVIIaIX - k11*TFVIIaIX;

ydot(13)=k6*TFVIIa*X - k17*TFVIIaX - k12*TFVIIaX;

ydot(14)=k6*VIIIaIXa*X - k18*VIIIaIXaX - k13*VIIIaIXaX;

ydot(15)=k9*VIIIaIXa - k7*VIIIa*IXa + k11*TFVIIaIX + k15*IX*Xa + k15*IX*VaXa;

ydot(16)=k10*VaXa - k8*Xa*Va + k12*TFVIIaX + k13*VIIIaIXaX;

ydot(17)=k10*VaXa - k8*Xa*Va + k1*V*Xa + k2*V*IIa + k2*V*mIIa;

ydot(18)=k9*VIIIaIXa - k7*VIIIa*IXa + k3*VIII*Xa + k4*VIII*IIa + k4*VIII*mIIa;

end 
