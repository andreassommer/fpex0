function [pl,ql,pr,qr] = FokkerPlanckBCFUN(xl,ul,xr,ur,t)
  % Boundary Condition Function for Fokker-Planck PDE
  %
  % PDE:  c(x,t,u,du/dx)*du/dt = x^{-m} d/dx [x^m * f(x,t,u,du/dx)] + s(x,t,u,du/dx)
  %
  % BCFUN:  p(x,t,u) + q(x,t)*f(x,t,u,du/dx) = 0
  %
  % INPUT:  xl --> left boundary coordinate (x_left)
  %         ul --> approximate solution at left boundary
  %         xr --> right boundary coordinate (x_right)
  %         ul --> approximate solution at right boundary
  %          t --> time
  %
  % OUTPUT: ql, pl --> p and q evaluated at x_left
  %         qr, pr --> p and q evaluated at x_right
  %
  %
  % Andreas Sommer, 2016-2022
  % andreas.sommer@iwr.uni-heidelberg.de
  % code@andreas-sommer.eu
  
  ql = 0;
  qr = 0;
  pl = ul;   % force u(x_left) to be zero
  pr = ur;   % force u(x_right) to be zero
  
end