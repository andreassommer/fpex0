% fpex0_cyber_logo.m
% Generiert ein technisch aufregendes "Cool Tech" Logo für fpex0
% Zeigt das De-Smearing von rechts (verschmiert) nach links (fpex0-Peak)

clear; clc; close all;

% 1. Quadratische Figure im Dark-Mode / Deep Space Look
% fig = figure('Color', [0.03 0.04 0.07], 'Position',[100, 100, 800, 800]);
backcolor = [0.03 0.04 0.10];
fig = figure('Color', backcolor, 'Position',[3200, -100, 1200, 1200]);
ax = axes('Position', [0.05, 0.05, 0.90, 0.90], 'Color', 'none');
hold(ax, 'on');

% 2. Mathematische Landschaft generieren (Optimiert für Rechtsdrall)
x = linspace(10, 100, 400);    % Temperatur-Analogon
beta = linspace(0, 15, 50);    % Heizrate (läuft nach hinten rechts)
[X, B] = meshgrid(x, beta);

% FPEX0-Modellierung: beta=0 ist extrem scharf, beta=15 ist breit und verschoben
mu = 35 + 2.2 * B;             % Drift: Peak wandert nach rechts hinten
sigma = 1.5 + 0.6 * B;         % Diffusion: Peak wird nach hinten breiter
amplitude = 120 ./ (sigma);    % Flächenerhalt: Vorne hoch, hinten flach

Z = amplitude .* exp(-((X - mu).^2) ./ (2 * sigma.^2));

% Ein bisschen künstliches Messrauschen NUR bei hohen Heizraten hinzufügen
noise = (rand(size(Z)) - 0.5) .* (B * 0.4);
Z = max(0, Z + noise);

% 3. Plotting mit "Glow-Effekt"
% Transluzente Oberfläche für den Tech-Look
h_surf = surf(ax, X, B, Z, 'EdgeColor', 'none', 'FaceColor', 'interp');
alpha(h_surf, 0.85);

% Markante Konturlinien (Wasserfall) oben drauf legen
h_mesh = mesh(ax, X, B, Z, 'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 1.2);

% Den ultimativen FPEX0-Zielpeak (Heizrate = 0) als leuchtendes Highlight fetten!
plot3(ax, x-0.3, zeros(size(x)), Z(1,:), 'Color', [1.0 0.0 0.4], 'LineWidth', 5);

% 4. Farb- und Lichtdesign ("Exciting Cyber")
% Benutzerdefinierte Map: Von Tiefblau/Violett (Rauschen) zu hellem Cyan/Neon-Grün (Peak)
custom_map = [
    linspace(0.1, 0.0, 100)', linspace(0.0, 0.4, 100)', linspace(0.3, 1.0, 100)';
    linspace(0.0, 0.0, 100)', linspace(0.4, 0.9, 100)', linspace(1.0, 0.9, 100)';
    linspace(0.0, 0.0, 100)', linspace(0.9, 1.0, 100)', linspace(0.9, 0.4, 100)'
];
% colormap(ax, custom_map);

% Dramatische Beleuchtung für plastische Kanten
view(ax, 304, 4.8);               % Blickwinkel: nach rechts hinten fluchtend
camlight('left');
lighting('gouraud');
material('shiny');

% 5. Text-Branding ("fpex0" im Hologramm-Style)
% Über dem scharfen Peak platziert
font1 = 'Berlin Sans FB';
font2 = 'Arial Rounded MT';
font3 = 'Cooper Black';
fontname = font3;
text(ax, -7, 13, max(Z(:))*0.85, 'FPEX_{0}', 'FontSize', 150, 'FontWeight', 'demi', ...
     'Color', [0.0 1.0 1.0], 'FontName', font3, 'HorizontalAlignment', 'left', 'FontAngle', 'normal');
text(ax, -7, 12.5, max(Z(:))*0.6, 'DSC', 'FontSize', 70, 'FontWeight', 'normal', ...
     'Color', [0.0 1.0 1.0], 'FontName', font3, 'HorizontalAlignment', 'left', 'FontAngle', 'normal');
text(ax, -7, 12.5, max(Z(:))*0.5, 'De-Smearing', 'FontSize', 70, 'FontWeight', 'normal', ...
     'Color', [0.0 1.0 1.0], 'FontName', font3, 'HorizontalAlignment', 'left', 'FontAngle', 'normal');


% 6. Clean-up für Icon-Export
axis(ax, 'tight');
zlim([-2 max(Z(:))*1.2]);
axis(ax, 'off');
set(fig, 'InvertHardcopy', 'off');

% 7. Als hochauflösendes Bild speichern
%exportgraphics(fig, 'fpex0_cyber_logo.png', 'Resolution', 300);
disp('Das "Exciting 3D-Cyber-Logo" wurde als fpex0_cyber_logo.png gespeichert!');
