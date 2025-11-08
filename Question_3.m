%% Clean up
clearvars; close all; clc;

%% Constants

contact_resistance = 10^-3;
T_initial = 20+273;
T_infinity = (180:10:220)+273;
eps = 0.8;
h = 35;
sigma = 5.67e-8; % Stefan-Boltzmann constant
k_p = 0.25;
k_m = 0.6;
k_b = 0.45;
c_p = 2000;
c_m = 3500;
c_b = 2700;
rho_p = 600;
rho_m = 1000;
rho_b = 1060;
T_final = 53+273;
dr = 0.3125/100;

%% Initial conditions
% rho = [rho_b,rho_m,rho_p];
% k = [k_b,k_m,k_p];
% Cp = [c_b,c_m,c_p]; 
dt=1;
t=0; 
% r = 0:dr:7.5/100; % radial distance from the center
r1=0:dr:0.0625;
r2=0.0625:dr:0.075; 
r = [r1,r2];
n = length(r); % number of spatial points
c = 1; %Iteration counter
% delta_T = 1;
% res = 1e-5;
Fo_p = ((k_p*dt)/(rho_p*c_p*dr^2));
Fo_m= ((k_m*dt)/(rho_m*c_m*dr^2));
Fo_b= ((k_b*dt)/(rho_b*c_b*dr^2));
% Initialize temperature array
T = ones(n,1) * T_initial;

for m = 1:numel(T_infinity)
    T_inf = T_infinity(m);

    % --- reset per ambient case ---
    T = ones(n,1)*T_initial;
    c = 1;
    time = 0;
    q_t = [];           % heat rate history (W)
    t_vec = []; 

while T(1,c) < T_final
%Outside boundary
i=26;
c=c+1;

T(i,c)=2*Fo_p*(((0.075 - (dr/2))/(0.075 - (dr/4)))*(T(i-1,c-1) - T(i,c-1))+((h*0.075*dr)/(k_p*(0.075 - (dr/4))))*(T_inf - T(i,c-1))+((sigma*eps*0.075*dr)/(k_p*(0.075-(dr/4))))*(T_inf^4 - T(i,c-1)^4)) + T(i,c-1);

% Pastry Interior Nodes

    for i= 23:25
        T(i,c) = T(i,c-1) + ( Fo_p * ( (r(i) - dr/2)*(T(i-1,c-1) - T(i,c-1)) + (r(i) + dr/2)*(T(i+1,c-1) - T(i,c-1)) ) / (r(i)) );
    end

%Pastry Side of Mushroom -> Pastry (Node 21)
T(22,c) = T(22,c-1)+((T(23,c-1)-T(22,c-1))*((2*Fo_p*(0.0625+(dr/2))/(0.0625+(dr/4))))) + (((2*0.0625*dt)/(contact_resistance*rho_p*c_p*dr*(0.0625+(dr/4))))*(T(21,c-1)-T(22,c-1)));

%Mushroom Side of Mushroom -> Pastry (Node 20)
T(21,c) = T(21,c-1)+((T(20,c-1)-T(21,c-1))*((2*Fo_m*(0.0625-(dr/2))/(0.0625-(dr/4))))) + (((2*0.0625*dt)/(contact_resistance*rho_m*c_m*dr*(0.0625-(dr/4))))*(T(22,c-1)-T(21,c-1)));

% Mushroom Interior Nodes
    for i= 18:20
        T(i,c) = T(i,c-1) + ( Fo_m * ( (r(i) - dr/2)*(T(i-1,c-1) - T(i,c-1)) + (r(i) + dr/2)*(T(i+1,c-1) - T(i,c-1)) ) / (r(i)) );
    end    

% Beef/mushroom Boundary
rho_avg = (rho_m+rho_b)/2;
c_avg = (c_m+c_b)/2;

T(17,c) = T(17,c-1) + (((k_b*dt*(0.05-(dr/2))*(T(16,c-1)-T(17,c-1))) + (k_m*dt*(0.05+(dr/2))*(T(18,c-1)-T(17,c-1))))/(rho_avg*c_avg*0.05*dr^2));

% Beef Interior Nodes
    for i=2:16 
        T(i,c) = T(i,c-1) + ( Fo_b * ( (r(i) - dr/2)*(T(i-1,c-1) - T(i,c-1)) + (r(i) + dr/2)*(T(i+1,c-1) - T(i,c-1)) ) / (r(i)) );
    end

% Center Node
T(1,c) = 4*Fo_b*T(2,c-1)+(1-(4*Fo_b))*T(1,c-1);

time = c*dt;
fprintf('%d',time)

q_rad = eps*sigma*2*pi*0.075*(T_inf^4-T(end,c)^4);
q_conv = h * 2 * pi * 0.075 * (T_inf - T(end,c));
q_t(c) = q_rad + q_conv;  

end


 % ---- Plot T vs r at selected time steps ----
figure(m); clf; hold on;

% valid indices (rounded and limited to >=1)
idx_now    = c;
idx_quarter = max(1, round(c/4));
idx_half    = max(1, round(c/2));
idx_3quarter = max(1, round(3*c/4));
idx_start   = 1;

plot(r, T(:,idx_start)   - 273.15, 'LineWidth', 2);
plot(r, T(:,idx_quarter) - 273.15, 'LineWidth', 2);
plot(r, T(:,idx_half)    - 273.15, 'LineWidth', 2);
plot(r, T(:,idx_3quarter)- 273.15, 'LineWidth', 2);
plot(r, T(:,idx_now)     - 273.15, 'LineWidth', 2);

xlabel('Radius (m)');
ylabel('Temperature (°C)');
title(sprintf('Temperature vs Radius at t = %.1f s', time));
legend('Start','25%','50%','75%','Current','Location','Best');
grid on;
drawnow;

q_total = sum(q_t(:));
fprintf('\nQ = %.6f J',q_total)
fprintf('\nT_S = %.2f C',T(end,c)-273)
end