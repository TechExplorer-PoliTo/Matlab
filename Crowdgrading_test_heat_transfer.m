%%
%STAINLESS STEEL 304
clear
clc
close all
format short

%Data
D = 2;  %m
L = 4;  %m
T0 = 1180;  %C
T0K = T0 +273.15; %K
Tair = 20;  %C
TairK = Tair+273.15; %K
ro_s= 7850;  %Kg/m^3
alfa_s = 3.6e-5;  %C^-1
cp_s = 502;  %J/(Kg K)
sigma_boltz = 5.670374419e-8; % W/(m^2*K^4)
F = 1;
eps = 0.8;
R = 8.3144598;
v_air = linspace(1,5,5);   %m/s


%%BIOT NUMBER in the worst case (T = T0, v = 5 m/s)
%1.1)Viscosity of air and thermal conductivity at 20C by using Chapman-Enskog, Eucken and Wilke's rule. Byron tables for Lennard-Jones parameters.

%%Chapman-Enskog, Eucken
%Oxygen
MO2 = 32; %g/mol = kg/mol
cp_O2 = 0.9189e3;  %J/KgK
sigmaO2 = 3.433;  %Angstrom
ljO2 = 113; %K   % Lennard jones
ljO2defq = TairK/ljO2; % target
ljO2def = [2.50 2.60];     
omegaO2 = [1.0933 1.0807];     % collision integral       
omegaO2_q = interp1(ljO2def, omegaO2, ljO2defq, 'linear');   % returns one value

mu_O2 = 2.6693*10^(-6)*(sqrt(MO2*TairK)/((sigmaO2)^2*omegaO2_q));%Pa*s
k_O2 = mu_O2*(cp_O2+ (5*R*10^3)/(4*MO2));


%Nitrogen
MN2 = 28.013; %g/mol
cp_N2 = 1.041e3;  %J/KgK
sigmaN2 = 3.667;  %Angstrom
ljN2 = 99.8; %K   % Lennard jones
ljN2defq = TairK/ljN2;         % target
ljN2def = [2.9 3];     
omegaN2 = [1.0482 1.0388];     % collision integral
omegaN2_q = interp1(ljN2def, omegaN2, ljN2defq, 'linear');  % returns one value

mu_N2 = 2.6693*10^(-6)*(sqrt(MN2*TairK)/((sigmaN2)^2*omegaN2_q));%Pa*s
k_N2= mu_N2*(cp_N2+ (5*R*10^3)/(4*MN2));                

%Argon
MAr = 39.948; %g/mol
cp_Ar = 0.520e3;  %J/KgK
sigmaAr = 3.432;  %Angstrom
ljAr = 122.4; %K   % Lennard jones
ljArdefq = TairK/ljAr;         % target
ljArdef = [2.3 2.4];     
omegaAr = [1.122 1.107];     % collision integral
omegaAr_q = interp1(ljArdef, omegaAr, ljArdefq, 'linear');  % returns one value

mu_Ar = 2.6693*10^(-6)*(sqrt(MAr*TairK)/((sigmaAr)^2*omegaAr_q));%Pa*s
k_Ar= mu_Ar*(cp_Ar+ (5*R*10^3)/(4*MAr));

%%Wilke's rule
M = [MO2 MN2 MAr];
mu = [mu_O2 mu_N2 mu_Ar];
k = [k_O2 k_N2 k_Ar];
y = [0.21 0.78 0.01];     %molar fractions of air

%phi
M_ratio = [];
for i = 1:3
    for j = 1:3
        row1(j) = M(i)/M(j);
        row2(j) = mu(i)/mu(j);
    end
   M_ratio(i,:) = row1;
   mu_ratio(i,:) = row2;
end

phi = 1/sqrt(8)*((1+M_ratio).^(-0.5)).*(1+(mu_ratio.^0.5).*((M_ratio').^0.25)).^2;

%viscosity of air 
for i = 1:3
    Num = y(i)*mu(i);
    for j = 1:3
        Den(j) = (y(j)*phi(i,j));
    end
    mu_cycle(i) = Num/sum(Den);
end

mu_mix = sum(mu_cycle); %Pa.s

%thermal conductivity of air
for i = 1:3
    Num = y(i)*k(i);
    for j = 1:3
        Den(j) = (y(j)*phi(i,j));
    end
    k_cycle(i) = Num/sum(Den);
end

k_mix = sum(k_cycle); %W/mK

%1.2) density of air applying the ideal gas law
M_mix = sum(y.*M)/1000;
ro = (101325*M_mix)/(R*TairK); %Kg/m^3

% 1.3) cp_mix
cp = [cp_O2 cp_N2 cp_Ar];
w = (y.*M)/(M_mix*1000);
cp_mix = sum(w.*cp);

%2)Dimensionless numbers (Reynolds, Prandtl)

%Reynolds
Re = (ro*D)/mu_mix*v_air;

%Prandtl
Pr =(mu_mix*cp_mix)/k_mix;

%%HEAT TRANSFER COEFFICIENT AND TEMPERATURE EVOLUTION 

dt = 5; %s

%Collision integral data

%oxygen 1453.15/113 = 12.8   ;   293.15/113 = 2.59
%Nitrogen 1453.15/99.8 = 14.56 ;   293.15/99.8 = 2.93
%Argon 1453.15/122.4 = 11.87   ;   293.15/122.4 = 2.39

% min = 2.39   max = 14.56

T_star_table = [2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90,3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00, 12.00, 14.00, 16.00];         
Omega_mu_table = [1.122, 1.107, 1.0933, 1.0807, 1.0691, 1.0583, 1.0482, 1.0388, 1.0300, 1.0217, 1.0139, 1.0066, 0.9996, 0.9931, 0.9868, 0.9809, 0.9753, 0.9699, 0.9647, 0.9598, 0.9551, 0.9506, 0.9462, 0.9420, 0.9380, 0.9341, 0.9304, 0.9268, 0.8962, 0.8727, 0.8538, 0.8380, 0.8244, 0.8018, 0.7836, 0.7683];

figCooling = figure('Name', 'Bloom Temperature Evolution', 'Color', 'w');
hold on

h_all_velocities = cell(1, length(v_air)); % Cell array to store h vectors of different lengths
figH = figure('Name', 'Heat Transfer Coefficient vs Temperature', 'Color', 'w');
hold on;


for i = 1:length(v_air)
    % Initialise time and temperature
    t = 0;
    T_current_K = T0K;

    % Arrays to store data
    time_array = t;
    T_array_K = T_current_K;
    Biot_evolution = [];
    h_evolution = [];

    %Time-Stepping Loop
    while T_current_K > (TairK+0.01)

        %%Chapman-Enksog
        % Oxygen 
        T_star_O2 = T_current_K / ljO2;
        omegaO2_T_q = interp1(T_star_table, Omega_mu_table, T_star_O2, 'linear', 'extrap'); 
        mu_O2_T = 2.6693e-6 * (sqrt(MO2 * T_current_K) / (sigmaO2^2 * omegaO2_T_q));

        % Nitrogen
        T_star_N2 = T_current_K / ljN2;
        omegaN2_T_q = interp1(T_star_table, Omega_mu_table, T_star_N2, 'linear', 'extrap');  
        mu_N2_T = 2.6693e-6 * (sqrt(MN2 * T_current_K) / (sigmaN2^2 * omegaN2_T_q));

        % Argon 
        T_star_Ar = T_current_K / ljAr;
        omegaAr_T_q = interp1(T_star_table, Omega_mu_table, T_star_Ar, 'linear', 'extrap');  
        mu_Ar_T = 2.6693e-6 * (sqrt(MAr * T_current_K) / (sigmaAr^2 * omegaAr_T_q));

        %%Wilke's rule 
        mu_T = [mu_O2_T mu_N2_T mu_Ar_T];

        for idx1 = 1:3
            for j = 1:3
                mu_ratio(idx1,j) = mu_T(idx1)/mu_T(j);
            end
        end

        phi2 = 1/sqrt(8) * ((1 + M_ratio).^(-0.5)) .* (1 + (mu_ratio.^0.5) .* (M_ratio'.^0.25)).^2;

        for idx1 = 1:3
            Num = y(idx1) * mu_T(idx1);
            for j = 1:3
                Den(j) = y(j) * phi2(idx1,j);
            end
            mu_cycle(idx1) = Num / sum(Den);
        end

        % Viscosity of air for this time step
        mu_T_mix = sum(mu_cycle);

        Nu_T = (0.4 * Re(i)^0.5 + 0.06 * Re(i)^(2/3)) * Pr^0.4 * (mu_mix / mu_T_mix)^0.25;
        h_T = (Nu_T * k_mix) / D;
        h_evolution(end+1) = h_T;

        %Specific Heat in J/(kg K)
        cp_cal_g_K = 0.1122 + 3.222e-5 * T_current_K;
        cp_current = cp_cal_g_K * 4184; 

        % Solid Conductivity
        k_solid_T = (8.116e-2 + 1.618e-4 * T_current_K) * 100;

        dT_diff = T_current_K - T0K;
        %alpha_V=3alpha_L  ( volumetric and linear)
        L_current = L * (1 + (alfa_s/3) * dT_diff);
        D_current = D * (1 + (alfa_s/3) * dT_diff);
        A = (2 * pi * (D_current/2)) * (L_current + D_current/2);
        % A = (2*pi*(D/2))*(L+D/2); %Surface Area of the cylinder
        % V = V0*(1+alfa_s*(T_current_K-T0K));
        V0 = pi*(D/2)^2*L;   %Volume of the cylinder
        V = V0 * (1 + alfa_s* dT_diff);
        Lc = V/A;     %Characteristic lenght of the cylinder

        %Biot  number calculation
        Bi_current = (h_T * Lc) / k_solid_T;
        Biot_evolution(end+1) = Bi_current;

        % Heat contribution
        q_conv = h_T * (TairK - T_current_K);
        q_rad = F * sigma_boltz * eps * (TairK^4 - T_current_K^4);

        dTdt = (A / (ro_s * cp_current * V0)) * (q_conv + q_rad);
        %dTdt = (A / (ro_s * cp_s * V0)) * (q_conv + q_rad); with constant
        % cp

        % Update time and temperature
        t = t + dt;
        T_current_K = T_current_K + dTdt * dt;

        % Store values
        time_array(end+1) = t;
        T_array_K(end+1) = T_current_K;
    end
    %Biot verification 
    fprintf('v = %.1f m/s | Max Bi: %.4f | Final Time: %.2f hours\n',v_air(i), max(Biot_evolution), t / 3600);

    h_all_velocities{i} = h_evolution;

    % Switch to the h figure and plot
    figure(figH); 
    plot(T_array_K(1:end-1) - 273.15, h_evolution, 'LineWidth', 1.5, 'DisplayName', sprintf('v = %.1f m/s', v_air(i)));
    
    %cooling
    figure(figCooling); 
    plot(time_array / 3600, T_array_K - 273.15, 'LineWidth', 1.5, 'DisplayName', sprintf('v = %.1f m/s', v_air(i)));

    % Final cooling time
    fprintf('Time required to cool at v = %.1f m/s: %.2f hours\n', v_air(i), t / 3600);
end

%Formatting for the Cooling time graph
figure(figCooling); 
title('Bloom Temperature Evolution');
xlabel('Time (hours)');
ylabel('Temperature (°C)');
yline(Tair, '--k', 'Ambient (20°C)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
legend('Location', 'northeast');
grid on;
hold off;

% Formatting for the Heat Transfer Coefficient graph
figure(figH);
title('Heat Transfer Coefficient (h) vs Temperature');
xlabel('Temperature (°C)');
ylabel('Heat Transfer Coefficient, h (W/m^2K)');
legend('Location', 'northeast');
grid on;
hold off;




