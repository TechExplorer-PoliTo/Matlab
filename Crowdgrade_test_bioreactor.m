clc
close all
clear 

%Data 
T_def= 37;  %°C
p_water = 993.36; %Kg m^-3
p_air = 1.138; %Kg m^-3
mu_w = 0.6922e-3;  %Pa*s
k = 0.03; %m^2/s^2
eps = 0.1; %m^2/s^3
mu_m = mu_w + (0.002*p_water*k^2)/eps;
T_interval = [30,40]; %°C
sigma_interval = [7.12*10^-2,6.96*10^-2];  %N/m
sigma_def = interp1(T_interval, sigma_interval,T_def, 'linear');
N = 100;
d_bubbles = logspace(log10(0.2e-3), log10(40e-3), N);   %mm
g = 9.8067;
T_interval_D = [20,40]; %°C
D_interval = [1.97*10^-9,3.24*10^-9];  %N/m
D_def = interp1(T_interval_D, D_interval,T_def, 'linear');


%%FIRST TASK
%Pure gas-liquid systems(absence of turbulence)

U_1 = zeros(1, N);
for i = 1:N
    % Initial guess from Stokes law
    U0 = (p_water - p_air) * g * d_bubbles(i)^2 / (18 * mu_w);
    U = U0;
    error_tol = 1; % Initialise with a high error 
    while error_tol >= 1e-6 
        Re = p_water * U * d_bubbles(i) / mu_w;
        Eo = (p_water - p_air) * g * d_bubbles(i)^2 / sigma_def;
        % C_D correlation 
        term1_1 = (16/Re) * (1 + 0.15 * Re^0.687);
        term1_2 =48/Re;
        term1 = min(term1_1,term1_2);
        term2 = (8/3) * (Eo / (Eo + 4));
        CD = max(term1, term2);
        % New U_T
        U_new = sqrt((4/3) * ((p_water - p_air)/p_water) * (g * d_bubbles(i) / CD));

        error_tol = abs(U_new - U) / U;
        U = U_new; 
    end
    U_1(i) = U;
end


%Contaminated gas-liquid systems(absence of turbulence)

U_2 = zeros(1, N);
for i = 1:N
    % Initial guess from Stokes law
    U0 = (p_water - p_air) * g * d_bubbles(i)^2 / (18 * mu_w);
    U = U0;
    error_tol = 1; % Initialise with a high error
    while error_tol >= 1e-6 
        Re = p_water * U * d_bubbles(i) / mu_w;
        Eo = (p_water - p_air) * g * d_bubbles(i)^2 / sigma_def;
        % C_D correlation
        term1 = (24/Re) * (1 + 0.15 * Re^0.687);
        term2 = (8/3) * (Eo / (Eo + 4));
        CD = max(term1, term2);
        % New U_T
        U_new = sqrt((4/3) * ((p_water - p_air)/p_water) * (g * d_bubbles(i) / CD));
        error_tol = abs(U_new - U) / U;
        U = U_new; 
    end
    U_2(i) = U;
end

%%SECOND TASK
%Pure gas-liquid systems(effect of turbulence)

U_1_turb = zeros(1, N);
for i = 1:N
    % Initial guess from Stokes law
    U0 = (p_water - p_air) * g * d_bubbles(i)^2 / (18 * mu_w);
    U = U0;
    error_tol = 1; % Initialise with a high error 
    while error_tol >= 1e-6 
        Re = p_water * U * d_bubbles(i) / mu_m;
        Eo = (p_water - p_air) * g * d_bubbles(i)^2 / sigma_def;
        % C_D correlation
        term1_1 = (16/Re) * (1 + 0.15 * Re^0.687);
        term1_2 =48/Re;
        term1 = min(term1_1,term1_2);
        term2 = (8/3) * (Eo / (Eo + 4));
        CD = max(term1, term2);
        % New U_T
        U_new = sqrt((4/3) * ((p_water - p_air)/p_water) * (g * d_bubbles(i) / CD));
        error_tol = abs(U_new - U) / U;
        U = U_new; 
    end
    U_1_turb(i) = U;
end


%Contaminated gas-liquid systems(effect of turbulence)

U_2_turb = zeros(1, N);
for i = 1:N
    % Initial guess from Stokes law
    U0 = (p_water - p_air) * g * d_bubbles(i)^2 / (18 * mu_w);
    U = U0;
    error_tol = 1; % Initialise with a high error
    while error_tol >= 1e-6 
        Re = p_water * U * d_bubbles(i) / mu_m;
        Eo = (p_water - p_air) * g * d_bubbles(i)^2 / sigma_def;
        % C_D correlation 
        term1 = (24/Re) * (1 + 0.15 * Re^0.687);
        term2 = (8/3) * (Eo / (Eo + 4));
        CD = max(term1, term2);
        % New U_T
        U_new = sqrt((4/3) * ((p_water - p_air)/p_water) * (g * d_bubbles(i) / CD));
        error_tol = abs(U_new - U) / U;
        U = U_new; 
    end
    U_2_turb(i) = U;
end


%%THIRD TASK
%Mass transfer coefficient by using the penetration theory (diffusion
%without reaction)
k_L_1      = 2 * sqrt((D_def .* U_1)      ./ (pi * d_bubbles));
k_L_2      = 2 * sqrt((D_def .* U_2)      ./ (pi * d_bubbles));
k_L_1_turb = 2 * sqrt((D_def .* U_1_turb) ./ (pi * d_bubbles));
k_L_2_turb = 2 * sqrt((D_def .* U_2_turb) ./ (pi * d_bubbles));


%%RESULTS

% Figure 1: Terminal Velocity
figure;
loglog(d_bubbles*10^3, U_1*10^2, 'b-', 'LineWidth', 1.5);        
hold on;
loglog(d_bubbles*10^3, U_2*10^2, 'r-', 'LineWidth', 1.5);       
loglog(d_bubbles*10^3, U_1_turb*10^2, 'g--', 'LineWidth', 1.5);   
loglog(d_bubbles*10^3, U_2_turb*10^2, 'm--', 'LineWidth', 1.5); 

title('Terminal Velocity of Bubbles');
xlabel('Bubble diameter [mm]');
ylabel('Terminal velocity [cm/s]');
legend('Pure System', 'Contaminated System', 'Pure System (Turbulent)', 'Contaminated System (Turbulent)', 'Location', 'southeast');
grid on;
hold off;

% Figure 2: Mass Transfer Coefficient
figure;
loglog(d_bubbles*10^3, k_L_1, 'b-', 'LineWidth', 1.5);
hold on;
loglog(d_bubbles*10^3, k_L_2, 'r-', 'LineWidth', 1.5);
loglog(d_bubbles*10^3, k_L_1_turb, 'g--', 'LineWidth', 1.5);
loglog(d_bubbles*10^3, k_L_2_turb, 'm--', 'LineWidth', 1.5);

title('Mass Transfer Coefficient (k_L) vs Bubble Diameter');
xlabel('Bubble diameter [mm]');
ylabel('Mass Transfer Coefficient, k_L [m/s]');
legend('Pure System', 'Contaminated System', 'Pure System (Turbulent)', 'Contaminated System (Turbulent)', 'Location', 'best');
grid on;
hold off;

