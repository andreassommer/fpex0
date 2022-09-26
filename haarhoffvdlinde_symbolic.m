function varargout = haarhoffvdlinde_symbolic(x,p)
   %               y = haarhoffvdlinde(x,p)
   %       [y, dydp] = haarhoffvdlinde(x,p)
   % [y, dydp, dydx] = haarhoffvdlinde(x,p)
   %
   % Haarhoff-Van der Linde function, also provides partial derivatives. SYMBOLIC VERSION
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
   %            dy --> derivate w.r.t. to p and x where
   %                   dy(1) = dyda0 --> derivate w.r.t. a0
   %                   dy(2) = dyda1 --> derivate w.r.t. a1
   %                   dy(3) = dyda2 --> derivate w.r.t. a2
   %                   dy(4) = dyda3 --> derivate w.r.t. a3
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


   persistent sym_dfdx sym_dfda0 sym_dfda1 sym_dfda2 sym_dfda3 sym_f
   persistent fun_dfdx fun_dfda0 fun_dfda1 fun_dfda2 fun_dfda3 fun_f
   
   % Generate derivatives ones
   if isempty(sym_dfdx) || (ischar(x) && strcmpi(x,'init'))
      fprintf('haarhoffvdlinde.m: Initializing partial derivatives... ');
      syms a0 a1 a2 a3 xx 
      sym_f = symfun(   ((a0 .* a2.^2) ./ (a1 .* a3 .* a2 .* sqrt(2*pi)) .* exp( - (xx-a1).^2 ./ (2 * a2.^2)) ) ...
                     ./ (1 ./ (exp(a1.*a3./a2.^2)-1) + 0.5 + 0.5*erf((xx-a1)./(a2*sqrt(2)))  )  ...
                     , [xx a0 a1 a2 a3]) ;  % order of variables
      sym_dfda0 = diff(sym_f,a0);
      sym_dfda1 = diff(sym_f,a1);
      sym_dfda2 = diff(sym_f,a2);
      sym_dfda3 = diff(sym_f,a3);
      sym_dfdx  = diff(sym_f,xx); % unused !
      fprintf('Transforming into matlab functions... ');
      fun_f     = matlabFunction(sym_f);
      fun_dfda0 = matlabFunction(sym_dfda0);
      fun_dfda1 = matlabFunction(sym_dfda1);
      fun_dfda2 = matlabFunction(sym_dfda2);
      fun_dfda3 = matlabFunction(sym_dfda3);
      fun_dfdx  = matlabFunction(sym_dfdx);
      fprintf('Done!\n');
   end

   
   % accessors
   a0 = p(1); a1 = p(2); a2 = p(3); a3 = p(4); 
   
   % first output argument: nominal values
   y = fun_f(x,a0,a1,a2,a3);
   varargout{1} = y;

   
   % derivatives requested?
   if (nargout >= 2)
      % derivatives w.r.t. parameters a0, a1, a2, a3
      dfda0 = fun_dfda0(x,a0,a1,a2,a3);
      dfda1 = fun_dfda1(x,a0,a1,a2,a3);
      dfda2 = fun_dfda2(x,a0,a1,a2,a3);
      dfda3 = fun_dfda3(x,a0,a1,a2,a3);
      
      % assemble output
      dfda0 = reshape(dfda0, [], 1);
      dfda1 = reshape(dfda1, [], 1);
      dfda2 = reshape(dfda2, [], 1);
      dfda3 = reshape(dfda3, [], 1);
      dydp = [dfda0 , dfda1 , dfda2 , dfda3];
      varargout{2} = dydp;
   end
   
   if (nargout >= 3)
      % derivative w.r.t. x
      dfdx = fun_dfdx(x,a0,a1,a2,a3);
      dfdx = reshape(dfdx, [], 1);
      varargout{3} = dfdx;
   end
   