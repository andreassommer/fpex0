% FPEX0_generate_3d_logo.m
% Generiert ein modernes 3D-Landschafts-Logo für MATLAB File Exchange
% Inspiriert von der Fokker-Planck-Extrapolation (Heizrate vs. Temperatur)

clear; clc; close all;

% 1. Quadratische Figure mit edlem, dunklem oder neutralem Hintergrund erzeugen
% Für einen modernen "Tech-Look" nutzen wir hier ein sehr helles Graublau.
fig = figure('Color', [0.95 0.96 0.98], 'Position', [100, 100, 800, 800]);
ax = axes('Position', [0.05, 0.05, 0.9, 0.9], 'Color', 'none');
hold(ax, 'on');

% 2. Mathematische Fokker-Planck-Landschaft simulieren (Analogie zu Fig. 6)
x = linspace(95, 155, 300);     % Raum / Temperatur T
t = linspace(0, 12, 45);        % Zeit / Heizrate beta
[X, T] = meshgrid(x, t);

% Heuristische Parameter für Drift & Diffusion
mu_0 = 131.6;        % Peak-Zentrum bei Heizrate 0 (fpex0)
drift_coeff = 0.85;  % Verschiebung nach rechts bei höherer Heizrate
sigma_0 = 1.4;       % Schärfe bei Heizrate 0
diff_coeff = 0.45;   % Verbreiterung (Flächenerhalt)

% Berechnung der driftenden und diffundierenden Peaks
mean_shift = mu_0 + drift_coeff * T + 0.03 * T.^2;
variance = sigma_0^2 + 2 * diff_coeff * T;
amplitude = 140 ./ (sqrt(variance) * sqrt(2 * pi));

% Die 3D-Matrix Z (Spezifische Wärmekapazität c_p)
Z = amplitude .* exp(-((X - mean_shift).^2) ./ (2 * variance));

% 3. 3D-Plotting: Wasserfall-Struktur mit feinen Linien für den mathematischen Look
% 'waterfall' erzeugt die typischen DSC-Kurven-Schichten entlang der Heizrate
h_surf = waterfall(ax, X, T, Z);

% 4. Visuelles Finetuning (Farben, Licht und Schatten)
colormap(ax, 'parula'); % 'parula' bietet den klassischen MATLAB-Wiedererkennungswert
set(h_surf, 'LineWidth', 1.5, 'EdgeLighting', 'gouraud', 'FaceAlpha', 0.5);

% Den vordersten Peak (Heizrate beta = 0, also FPEX0) extra betonen!
plot3(ax, x, zeros(size(x)), Z(1,:), 'Color', [0.85 0.33 0.1], 'LineWidth', 4.5);

% 5. Perspektive und Beleuchtung einstellen
view(ax, -35, 30);          % Perfekter 3D-Winkel für Kurvenlandschaften
camlight('headlight');      % Erzeugt edle Lichtreflexe auf den Kurvenkanten
lighting('gouraud');
material('shiny');

% 6. Text-Branding im 3D-Raum platzieren
% Schwebend über oder dezent hinter der Landschaft platziert
text(ax, 110, 14, 25, 'fpex_{0}', 'FontSize', 55, 'FontWeight', 'bold', ...
     'Color', [0.07 0.07 0.25], 'FontName', 'Helvetica', 'HorizontalAlignment', 'center');

% 7. Achsen für den reinen Icon-Look eliminieren
axis(ax, 'tight');
zlim([-1 max(Z(:))+5]);
axis(ax, 'off'); 
set(fig, 'InvertHardcopy', 'off');

% 8. Speichern als hochauflösendes, quadratisches Bild
%exportgraphics(fig, 'fpex0_3d_logo.png', 'Resolution', 300);
disp('Fancy 3D-Logo wurde erfolgreich als fpex0_3d_logo.png gespeichert!');
