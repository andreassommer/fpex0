function [c,f,s] = FokkerPlanckPDEFUN_parametrized(x,t,u,dudx,p)
  % Fokker-Planck PDE function, with parametrized drift/diffusion/drain
  %
  % PDE:  c(x,t,u,du/dx)*du/dt = x^{-m} d/dx [x^m * f(x,t,u,du/dx)] + s(x,t,u,du/dx)
  %
  % INPUT:  x --> current position
  %         t --> current time
  %         u --> approximate solution
  %     du/dx --> approximate u'
  %         p --> cell array of parameter functions {drift(t), diffusion(t), drain(t)}
  %
  % OUTPUT: c --> coefficent function c as above
  %         f --> coefficent function f as above (flux)
  %         s --> coefficent function s as above (source)
  %
  % Andreas Sommer, 2016-2022
  % andreas.sommer@iwr.uni-heidelberg.de
  % code@andreas-sommer.eu

  mu    = p{1};  % drift
  D     = p{2};  % diffusion
  drain = p{3};  % drain
  
  % coefficient of lhs is 1
  c = 1;
  
  % flux is defined by Fokker-Planck
  f = - mu(t) * u  +  D(t) * dudx;
  
  % source term --> small drain!
  if length(p)==3
     %s = drain(t)* max(0,u);         % linear drain 
     s = - (exp(drain(t) * u) - 1) ;  % exponential drain
  else
     s = 0;                           % nodrain
  end
  
end

