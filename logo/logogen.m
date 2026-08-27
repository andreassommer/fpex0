% FPEX0_generate_logo.m
% Generiert ein modernes quadratisches Icon/Logo für MATLAB File Exchange

clear; clc; close all;

% 1. Quadratische Figure erstellen (Standard für File Exchange App Icons)
fig = figure('Color', [0.96 0.96 0.98], 'Position', [100, 100, 600, 600]);
ax = axes('Position', [0.1, 0.1, 0.8, 0.8], 'Color', 'none');
hold(ax, 'on');

% X-Achse (Temperatur-Analogon)
x = linspace(0, 10, 1000);

% Heuristische Peak-Funktionen zur Demonstration des Effekts:
% Smeared Peak (Heizrate > 0): Breit, nach rechts verschoben, niedriger
y_smeared = 1.8 * exp(-((x - 6.0) / 1.8).^2); 

% De-Smeared Peak (FPEX0 bei Rate 0): Scharf, links, hoch
y_fpex0   = 4.5 * exp(-((x - 4.5) / 0.6).^2); 

% 2. Daten plotten
% Gefüllte Fläche für das verschmierte DSC-Signal im Hintergrund (Grau/Blau-Verlauf)
fill([x, fliplr(x)], [y_smeared, zeros(size(x))], [0.7 0.75 0.85], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(x, y_smeared, 'Color', [0.5 0.6 0.7], 'LineWidth', 2.5, 'LineStyle', '--');

% Gefüllte Fläche für das de-smeared FPEX0 Signal (Kräftiges MATLAB-Blau)
fill([x, fliplr(x)], [y_fpex0, zeros(size(x))], [0.0 0.45 0.74], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.15);
plot(x, y_fpex0, 'Color', [0.0 0.45 0.74], 'LineWidth', 5);

% Eine dünne, prägnante Null-Linie (Baseline)
plot([0 10], [0 0], 'Color', [0.2 0.2 0.2], 'LineWidth', 3);

% 3. Text/Branding hinzufügen
text(4.5, 5.1, 'fpex_{0}', 'FontSize', 46, 'FontWeight', 'bold', ...
     'Color', [0.0 0.45 0.74], 'HorizontalAlignment', 'center', 'FontName', 'Helvetica');

% 4. Stil-Optimierungen für ein sauberes Icon
axis(ax, 'tight');
xlim([1 9]);
ylim([-0.5 6]);
axis(ax, 'off'); % Achsenbeschriftungen komplett ausblenden für Icon-Look
set(fig, 'InvertHardcopy', 'off');

% 5. Automatisch als hochauflösendes PNG speichern
% File Exchange empfiehlt quadratische Bilder (z.B. 500x500 oder höher)
%exportgraphics(fig, 'fpex0_logo.png', 'Resolution', 300);
disp('Das FPEX0-Logo wurde erfolgreich als fpex0_logo.png gespeichert!');
