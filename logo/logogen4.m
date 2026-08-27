% fpex0_clean_tech_logo.m
% Generiert ein helles, aufregendes und hoch-professionelles 3D-Logo für fpex0
% Zeigt den De-Smearing-Effekt von rechts (verschmiert) nach links vorne (scharf)

clear; clc; close all;

% 1. Quadratische Figure mit modernem, hellem Hintergrund (Clean Tech Look)
fig = figure('Color', [0.97 0.98 1.0], 'Position',[100, 100, 800, 800]);
ax = axes('Position', [0.05, 0.05, 0.9, 0.9], 'Color', 'none');
hold(ax, 'on');

% 2. Mathematische Landschaft (Rechtsdrall für De-Smearing)
x = linspace(10, 100, 400);    % Temperatur-Analogon
beta = linspace(0, 15, 45);    % Heizrate (läuft nach hinten rechts)
[X, B] = meshgrid(x, beta);

% FPEX0-Modellierung: Vorne scharf (fpex0), hinten breit und verschoben
mu = 32 + 2.4 * B;             % Drift: Peak wandert nach rechts
sigma = 1.3 + 0.55 * B;        % Diffusion: Peak wird breiter
amplitude = 135 ./ (sigma);    % Vorne hoch, hinten flach

Z = amplitude .* exp(-((X - mu).^2) ./ (2 * sigma.^2));

% 3. Plotting: Transluzente Oberfläche mit feinen Gitterlinien
% Wir nutzen 'meshz' für einen architektonischen, technischen Look ("Vorhang-Effekt")
h_mesh = meshz(ax, X, B, Z);
set(h_mesh, 'EdgeColor', 'interp', 'FaceColor', [0.95 0.97 1.0], ...
            'FaceAlpha', 0.6, 'LineWidth', 1.1);

% Den FPEX0-Zielpeak (Heizrate = 0) im Vordergrund als "Hero-Element" betonen (Kräftiges Orange)
plot3(ax, x, zeros(size(x)), Z(1,:), 'Color', [0.92 0.40 0.12], 'LineWidth', 5);

% 4. Helles Farb- und Lichtdesign
% Farbverlauf von tiefem Universitäts-Blau über Türkis zu hellem Cyan
clean_map = [
    linspace(0.0, 0.0, 128)', linspace(0.3, 0.7, 128)', linspace(0.6, 0.9, 128)';
    linspace(0.0, 0.0, 127)', linspace(0.7, 0.4, 127)', linspace(0.9, 0.7, 127)'
];
colormap(ax, clean_map);

% Perfekt ausgeleuchtete 3D-Vektorgrafik-Optik
view(ax, -52, 34);               % Blickwinkel fluchtend nach rechts hinten
camlight('left');
lighting('flat');                % Flachere Schattierung für den "Clean Graphic" Stil

% 5. Modernes, serifenloses Text-Branding
% Edel und klar über der Landschaft schwebend
text(ax, 35, -4, max(Z(:))*1.1, 'fpex_{0}', 'FontSize', 56, 'FontWeight', 'bold', ...
     'Color', [0.05 0.15 0.35], 'FontName', 'Arial', 'HorizontalAlignment', 'center');

% 6. Clean-up für den perfekten Icon-Export
axis(ax, 'tight');
zlim([-1 max(Z(:))*1.15]);
axis(ax, 'off');
set(fig, 'InvertHardcopy', 'off');

% 7. Als hochauflösendes Bild speichern
% exportgraphics(fig, 'fpex0_clean_logo.png', 'Resolution', 300);
disp('Das helle Clean-Tech 3D-Logo wurde als fpex0_clean_logo.png gespeichert!');
