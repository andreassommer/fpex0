function varargout = haarhoffvdlinde(x,p)
   %               y = haarhoffvdlinde(x,p)
   %       [y, dydp] = haarhoffvdlinde(x,p)
   % [y, dydp, dydx] = haarhoffvdlinde(x,p)
   %
   % Haarhoff-Van der Linde function, also provides partial derivatives.
   %
   % See 
   %  Haarhoff, Van der Linde 1966: Concentration dependence of Elution Curves in Non-Ideal Gas Chromatography
   %  Le Saux, Varenne, Gareil 2005: Peak shape modeling by HVL function for the determination of correct migrationtimes
   %
   % INPUT:      x --> evaluation points
   %             p --> parameter vector where
   %                   p(1) = a0 --> peak area
   %                   p(2) = a1 --> center of Gaussian component
   %                   p(3) = a2 --> standard deviation of the Gaussian component
   %                   p(4) = a3 --> peak distortion
   %
   % OUTPUT:     y --> y(x) values at position x
   %            dy --> derivative w.r.t. to p where
   %                   dy(1) = dyda0 --> derivate w.r.t. a0
   %                   dy(2) = dyda1 --> derivate w.r.t. a1
   %                   dy(3) = dyda2 --> derivate w.r.t. a2
   %                   dy(4) = dyda3 --> derivate w.r.t. a3
   %            dx --> derivative w.r.t. to x
   %
   % CONSTRAINTS:  a0 > 0,  a2 > 0
   %               No check is done on constraints!
   %
   % Andreas Sommer, Dec2021
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu
   %
   %
   % Copyright 2016-2022, Andreas Sommer  code@andreas-sommer.eu
   %
   % Copying and distribution of this file, with or without modification, are permitted in any medium without royalty,
   % provided the copyright notice and this notice are preserved. This file is offered as-is, without any warranty.

   
   % accessors
   a0 = p(1); a1 = p(2); a2 = p(3); a3 = p(4); 
   
   % reshape x to single column, store original dimension
   x_shape = size(x);
   x = reshape(x, [], 1);
   
   % common subexpressions (used in f and df)
   sqrt2     = sqrt(2);
   sqrt2pi   = sqrt(2*pi);
   a1mX      = a1 - x;
   a1mX2     = a1mX .^ 2;
   exp_a1mX2 = exp(-a1mX2);
   A         = 1 ./ ( 1/(exp((a1*a3)/a2^2)-1) + 0.5 - erf(a1mX/(a2*sqrt2))/2 );  % bad scaling and cancellation!
   % For higher accuracy, one might think using erf or erfc=1-erf whenever appropriate
   B         = exp_a1mX2 .^ (1/(2*a2^2));
   AB        = A .* B ;
   a0AB      = a0 * AB;
   a0a2AB    = a2 * a0AB;
   a1a3sqrt2pi = a1*a3*sqrt2pi;
   
   
   % FIRST OUTPUT: nominal values
   y = a0a2AB / a1a3sqrt2pi;
   y = reshape(y, x_shape);    % reshape to x input shape
   varargout{1} = y;

   
   % derivatives requested?
   if (nargout > 1)
      
      % common subexpressions (used only by df)
      Asq   = A.^2;
      C     = exp_a1mX2 .^ (1/(a2^2*sqrt2^2));
      D     = exp((a1*a3)/a2^2);
      AsqB  = Asq .* B;
      Dm1sq = (D - 1).^2 ;
      D_by_Dm1sq = D ./ Dm1sq;
      a0AsqB     = a0 * AsqB;
      a0a2AsqB   = a2 * a0AsqB;  
      
      % derivatives w.r.t. parameters a0, a1, a2, a3
      dfda0 = AB * a2 / a1a3sqrt2pi;
      dfda1 = (a0a2AsqB .* ( a3/a2^2 * D_by_Dm1sq  + C ./ (a2*sqrt2*pi^0.5)) ) / a1a3sqrt2pi            ...
              - a0AB  .* ( 2*a1mX / (2*a2*a1a3sqrt2pi) - a2 / (a1*a1a3sqrt2pi) );
      dfda2 = a0AB / a1a3sqrt2pi  +  (a0*AB.*a1mX2) / (a2^2*a1a3sqrt2pi)                              ...
              - (a0a2AsqB .* ( (2*a1*a3*D)./(a2^3*Dm1sq) + (C .* a1mX)/(a2^2*sqrt2*pi^0.5))) / a1a3sqrt2pi;
      dfda3 = a0AsqB .* D_by_Dm1sq ./ (a2*a3*sqrt2pi)  -  a0a2AB/(a3*a1a3sqrt2pi);
      
      % assemble output
      dfda0 = reshape(dfda0, [], 1);
      dfda1 = reshape(dfda1, [], 1);
      dfda2 = reshape(dfda2, [], 1);
      dfda3 = reshape(dfda3, [], 1);
      dydp = [dfda0 , dfda1 , dfda2 , dfda3];
      
      varargout{2} = dydp;
   end
   
   
   % derivatives w.r.t. x requested?
   if (nargout >= 3)
      
      % derivates w.r.t. x
      dfdx = (2*a0AB .* a1mX) / (2*a2*a1a3sqrt2pi)  -  (AsqB.*C*a0) / (a1a3sqrt2pi*sqrt2*pi^0.5);      
      dfdx = reshape(dfdx, [], 1);
      
      varargout{3} = dfdx;
   end

   
   % finito
   return
   
end



% Code to generate some common subexpressions
%  syms a0 a1 a2 a3 xx sqrt2pi sqrt2
%  sym_f =  ((a0 .* a2.^2) ./ (a1 .* a3 .* a2 .* sqrt2pi) .* exp( - (xx-a1).^2 ./ (2 * a2.^2)) ) ...
%           ./ (1 ./ (exp(a1.*a3./a2.^2)-1) + 0.5 + 0.5*erf((xx-a1)./(a2*sqrt2))  ) ;
%  sym_dfda0 = diff(sym_f,a0);
%  sym_dfda1 = diff(sym_f,a1);
%  sym_dfda2 = diff(sym_f,a2);
%  sym_dfda3 = diff(sym_f,a3);
%  sym_dfdx  = diff(sym_f,xx);
%  
%  ss0 = [sym_f ; sym_dfda0 ; sym_dfda1 ; sym_dfda2 ; sym_dfda3 ; sym_dfdx ];
%  
%  [ss1, A] = subexpr(ss0, 'A');
%  [ss2, B] = subexpr(ss1, 'B');
%  [ss3, C] = subexpr(ss2, 'C');
%  [ss4, D] = subexpr(ss3, 'D');
 