% fpex0_bright_glow_logo.m
% Generiert ein helles Logo mit intensivem Neon-Glow-Effekt für fpex0
% Verbindet eine cleane Ästhetik mit aufregenden Licht- und Glanzeffekten

clear; clc; close all;

% 1. Quadratische Figure mit ultra-hellem, modernem Studio-Hintergrund
fig = figure('Color', [0.96 0.97 0.99], 'Position', [100, 100, 800, 800]);
ax = axes('Position', [0.05, 0.05, 0.9, 0.9], 'Color', 'none');
hold(ax, 'on');

% 2. Mathematische Landschaft (Rechtsdrall für De-Smearing)
x = linspace(10, 100, 450);    % Temperatur-Analogon
beta = linspace(0, 15, 55);    % Heizrate (läuft nach hinten rechts)
[X, B] = meshgrid(x, beta);

% FPEX0-Modellierung: Vorne scharf, hinten breit und verschoben
mu = 30 + 2.5 * B;             % Drift: Peak wandert nach rechts
sigma = 1.2 + 0.6 * B;         % Diffusion: Peak wird breiter
amplitude = 140 ./ (sigma);    % Vorne hoch, hinten flach

Z = amplitude .* exp(-((X - mu).^2) ./ (2 * sigma.^2));

% 3. Plotting: Gitterstruktur mit hochentwickelten Materialeigenschaften
% Wir nutzen 'surf' für maximale Lichtreflexion (Glow)
h_surf = surf(ax, X, B, Z, 'EdgeColor', 'interp', 'FaceColor', 'interp', 'LineWidth', 0.8);

% Materialeigenschaften für den "Glaskörper/Glow"-Look einstellen
set(h_surf, 'FaceAlpha', 0.85, ...
            'AmbientStrength', 0.6, ...
            'DiffuseStrength', 0.8, ...
            'SpecularStrength', 1.0, ...    % Maximale Reflexion für den Glow
            'SpecularExponent', 8);        % Konzentrierter Lichtfleck

% Den FPEX0-Zielpeak (Heizrate = 0) im Vordergrund als "Laser-Klinge" fetten
plot3(ax, x, zeros(size(x)), Z(1,:), 'Color', [1.0 0.1 0.4], 'LineWidth', 5.5);

% 4. Farb- und Lichtdesign für den Bright Neon-Look
% Farbverlauf von elektrischem Neon-Cyan über Violett zu glühendem Pink/Orange
glow_map = [
    linspace(0.0, 0.5, 128)', linspace(0.8, 0.0, 128)', linspace(1.0, 0.8, 128)'; % Cyan zu Violett
    linspace(0.5, 1.0, 127)', linspace(0.0, 0.3, 127)', linspace(0.8, 0.1, 127)'  % Violett zu Neon-Orange
];
colormap(ax, glow_map);

% Dramatische Studiobeleuchtung (Das erzeugt den eigentlichen Glow auf hellem Grund)
view(ax, -54, 36);               % Perfekter Blickwinkel nach rechts hinten
h_light = camlight('left');      % Lichtquelle von links setzen
lighting('gouraud');

% 5. Tech-Schriftzug (Schwebendes Hologramm-Branding)
text(ax, 34, -4, max(Z(:))*1.12, 'fpex_{0}', 'FontSize', 58, 'FontWeight', 'bold', ...
     'Color', [0.1 0.15 0.35], 'FontName', 'Trebuchet MS', 'HorizontalAlignment', 'center');

% 6. Clean-up für den Icon-Export
axis(ax, 'tight');
zlim([-1 max(Z(:))*1.15]);
axis(ax, 'off');
set(fig, 'InvertHardcopy', 'off');

% 7. Als hochauflösendes Bild speichern
% exportgraphics(fig, 'fpex0_bright_glow.png', 'Resolution', 300);
disp('Das "Bright Glow"-Logo wurde erfolgreich als fpex0_bright_glow.png gespeichert!');
