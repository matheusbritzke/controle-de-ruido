% PLACA PERFURADA
clear; close all; clc;

function [Zc, kc]=DB(f,rho0,c0,sigma)
% Impedância característica do material
Zc = (rho0 * c0).* ( (1 + 9.08 .* ((1e3.*f)./sigma).^(-0.75)) - 1j * 11.9 .* ((1e3.*f)./sigma).^(-0.73) );

k0 = (2*pi.*f)./c0;
% Número de onda característico:
kc = k0.* ( (1 + 10.8.*((1e3.*f)./sigma).^(-0.70)) - 1j * (10.3 .* ((1e3.*f)./sigma).^(-0.59)) );
end

% FREQUENCIAS 
f_i = 20; % [Hz]Frequencia inicial
f_f = 10000; % [Hz] Frequencia final

f_s = (f_f - f_i);
f = linspace(f_i,f_f,f_s); % [Hz] Vetor de frequencia linear
omega = 2*pi*f; % [rad/s] Vetor de frequencia angular

% PROPRIEDADES AR
rho0 = 1.21; % [kg/m^3] densidade do ar
c0 = 343; % [m/s] velocidade do som no ar
k0 = omega/c0; % [rad/m] numero de onda
viscos = 1.84e-5; % [Pa.s] Viscosidade do ar (usado em denssup)


%PROPRIEDADES DO ABSORVEDOR
%POROSO
sigma = 20000; % [Ns/m^4] Resistividade ao fluxo do material poroso
d = 6e-2; % [m] Espessura do material poroso

D = .12; % [m] comprimento total da cavidade do absorvedor
%PLACA
a = 0.0025; % [m] raio do furo da placa
b = 0.014; % [m] espaçamento entre os furos 
l = 0.001; % [m] espessura da placa perfurada 
lc = l + 1.7*a; % [m] comprimento corrigido (impedância de radiação)
psi = (pi*(a*2)^2)/(4*b^2); % razão de perfuração da placa
denssup = rho0.*(1/psi).*(lc+(((8*viscos./omega).*(1+(l/2*a))).^0.5)); % densidade superficial


[Zc, kc] = DB(f, rho0, c0, sigma);
Zpp = 1j.*omega.*denssup;
Zsar =  (-1j * rho0 * c0) ./ (tan(k0 .* (D-d)));
Zsi = ((-1j .* Zc .* Zsar ./ tan(kc.*d)) + Zc.^2) ./ (Zsar - ((1j .* Zc) ./ tan(kc.*d)));
Zs = Zpp + Zsi;

% Coeficiente de absorção alpha                   
reflex = (Zs - (rho0*c0)) ./ (Zs + (rho0*c0)); % Coeficiente de reflexão
alpha = 1 - (abs(reflex).^2);


semilogx(f, alpha, 'LineWidth', 2);
title('Placa Perfurada','FontSize', 13');
%subtitle(sprintf('\\sigma = %.0f N·s/m^4 | d_por = %.3f m | d_mem = %.4f m' , sigma, d, dm),'FontSize', 11)
% Configurações finais do gráfico
set(gca, 'XScale', 'log')
set(gca, 'TickLabelInterpreter', 'tex') 
xlabel('Frequência (Hz)')
ylabel('\alpha (Coeficiente de Absorção)')
grid on
axis([f(1), f(end), 0, 1.05])