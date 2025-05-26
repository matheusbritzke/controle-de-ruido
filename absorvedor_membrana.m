% ABSORVEDOR DE MEMBRANA
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

f_s = (f_f - f_i)*2;
f = linspace(f_i,f_f,f_s); % [Hz] Vetor de frequencia linear
omega = 2*pi*f; % [rad/s] Vetor de frequencia angular

% PROPRIEDADES AR
rho0 = 1.21; % [kg/m^3] densidade do ar
c0 = 343; % [m/s] velocidade do som no ar
k0 = omega/c0; % [rad/m] numero de onda

%PROPRIEDADES DO ABSORVEDOR

sigma = 25000; % [Ns/m^4] Resistividade ao fluxo do material poroso
D = .20; % [m] comprimento total da cavidade do absorvedor
d = 7e-2; % [m] Espessura do material poroso
dm   = 2e-4; % [m] Espessura da membrana
rhom = 1300; % [kg/m^2] Densidade do material da membrana
denssup = dm*rhom; % [kg/m^2] Densidade superficial da membrana

Zm  = 1j .* omega .* denssup; % Impedância da membrana
[Zp, kp] = DB(f, rho0, c0, sigma); % Impedância característica do material poroso
Zsp = Zp .* ((cosh(kp.*d))./(sinh(kp.*d))); % Impedância no topo da camada de material poroso         
Zsi1 = (-1j.*Zsp.*rho0.*c0.*(1./tanh(k0.*(D-d)))) + ((rho0.*c0).^2);  % Impedância no topo da camada de ar
Zsi2 = Zsp - (1j.*rho0.*c0.*(1./tanh(k0.*(D-d))));

Zsi= Zsi1 ./ Zsi2 ;

Zs = Zm + Zsi; % Impedância de superfície do absorvedor

% Coeficiente de absorção alpha                   
reflex = (Zs - (rho0*c0)) ./ (Zs + (rho0*c0)); % Coeficiente de reflexão
alpha = 1 - (abs(reflex).^2);


semilogx(f, alpha, 'LineWidth', 2);
title('Absorvedor de Membrana','FontSize', 13');
subtitle(sprintf('\\sigma = %.0f N·s/m^4 | d_por = %.3f m | d_mem = %.4f m' , sigma, d, dm),'FontSize', 11)
% Configurações finais do gráfico
set(gca, 'XScale', 'log')
set(gca, 'TickLabelInterpreter', 'tex') 
xlabel('Frequência (Hz)')
ylabel('\alpha (Coeficiente de Absorção)')
grid on
axis([f(1), f(end), 0, 1.05])