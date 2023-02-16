%% This script contains a mouving boundary diffusion problem using an analytical solution

clear variables;

%% INPUT

C_0 = 0.3; %initial concentration
C_1 = 0.7; %boundary concentration
t = 60000*60*60; %time in s
dx = 1; %spatial resolution in um
l = 1000; %profile length in um
D = 1E-4; %diffusion coefficient in um²/s
v = 1E-6; %veolocity of boundary in um/s

%% CREATE DISTANCE ARRAY

%calculate number of nodes & allocate array
nodes = l/dx;
x = zeros(nodes,1);
term1 = zeros(nodes,1);
term2 = zeros(nodes,1);
term3 = zeros(nodes,1);
C = zeros(nodes,1);
x(1) = 0;

C_initial(1) = C_1;
C_initial(2:nodes) = C_0;
C_initial = C_initial';

for i = 2:nodes
    x(i) = x(i-1)+dx;
end

%% ANALYTICAL SOLUTION

%TERM I, II & III
for i = 1:nodes
    term1(i) = erfc((x(i)+v*t)/(2*sqrt(D*t)));
    term2(i) = exp((-v*x(i))/D);
    term3(i) = erfc((x(i)-v*t)/(2*sqrt(D*t)));
    C(i) = C_0+0.5*(C_1-C_0)*(term1(i)+term2(i)*term3(i));
end

%% ADDITION TO SHOW PROGRESS OF BOUNDARY

dist_from0_boundary = v*t;
x_moved = x+dist_from0_boundary;
x_add = (0:dx:(x_moved(1)-dx));
x_add = x_add';
x_extended = [x_add;x_moved];
C_add(1:numel(x_add)) = C_1;
C_add = C_add';
C_extended = [C_add;C];
C_initial_extended(1:numel(x_extended)) = C_0;
C_initial_extended(1) = C_1;
C_initial_extended = C_initial_extended';


%% PLOT RESULTS

figure
subplot(2,1,1)
plot(x,C_initial,'--','LineWidth',3,'Color','b');
hold on;
plot(x,C,'-','LineWidth',3,'Color','r');
ylim([0 1]);
legend('initial','calculation')
subplot(2,1,2)
plot(x_extended,C_initial_extended,'--','LineWidth',3,'Color','b')
hold on;
plot(x_extended,C_extended,'-','LineWidth',3,'Color','r')
ylim([0 1]);
legend('initial','calculation')
hold off;