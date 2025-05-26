% CLAUSURA
clear; close all; clc;

% FREQUENCIAS 
f_i = 20; % [Hz]Frequencia inicial
f_f = 10000; % [Hz] Frequencia final

f_s = (f_f - f_i);
f = linspace(f_i,f_f,f_s); % [Hz] Vetor de frequencia linear
omega = 2*pi*f; % [rad/s] Vetor de frequencia angular

% AR
rho0 = 1.21; 
c0 = 343;

% CLAUSURA
Li_xyz = [.5,.4,.3]; %[m] dimensões internas

% MATERIAL (MDF)
d_m = 0.02; % [m] espessura do material da clausura
rho_m = 750; % [kg/m^3] densidade do material
denssup = rho_m * d_m; % [kg/m^2] densidade superficial
alhpa_m = 0.1; % coeficiente de absorção (varia com a frequência mas usarei uma média)

% Perda de transmissão
TL = 10 * log10(1 + ( (omega .* denssup)./(2*rho0*c0) ) );

% Áreas da clausura
Si = 2*(Li_xyz(1)*Li_xyz(2) + Li_xyz(1)*Li_xyz(3) + Li_xyz(2)*Li_xyz(3)); % Area interna
Le_xyz = Li_xyz+ 2*d_m; % Dimensões externas
Se = 2*(Le_xyz(1)*Le_xyz(2) + Le_xyz(1)*Le_xyz(3) + Le_xyz(2)*Le_xyz(3)); % Area externa

% Fator de coreção
C = 10 * log10 ( ( 0.3 + Se * (1-alhpa_m)/(Si*alhpa_m) ) );
% se alpha fosse dependente da frequencia, C também seria.

IL = TL - C;



figure; 
% semilogx(f, IL); 
plot(f, IL, LineWidth=2);
title('Perda por inserção'); 
xlabel('Frequência [Hz]');
ylabel('Perda por inserção IL [dB]'); 
grid on;             
grid minor;
xlim([f(1), f(end)])