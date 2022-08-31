function [c,f,s] = FokkerPlanckPDEFUN(x,t,u,dudx,p)
  % Fokker-Planck PDE function
  %
  % PDE:  c(x,t,u,du/dx)*du/dt = x^{-m} d/dx [x^m * f(x,t,u,du/dx)] + s(x,t,u,du/dx)
  %
  % INPUT:  x --> current position
  %         t --> current time
  %         u --> approximate solution
  %     du/dx --> approximate u'
  %         p --> parameter vector
  %
  % OUTPUT: c --> coefficent function c as above
  %         f --> coefficent function f as above (flux)
  %         s --> coefficent function s as above (source)
  %
  % Andreas Sommer, 2016-2022
  % andreas.sommer@iwr.uni-heidelberg.de
  % code@andreas-sommer.eu

  mu = p(1);  % drift
  D = p(2);   % diffusion
  
  % coefficient of lhs is 1
  c = 1;
  
  % flux is defined by Fokker-Planck
  f = - mu * u  +  D * dudx;
  
  % source term --> small drain!
  if length(p)==3
     %s = p(3)* max(0,u);         % linear drain 
     s = - (exp(p(3) * u) - 1) ;  % exponential drain
  else
     s = 0;                       % nodrain
  end
  
end

