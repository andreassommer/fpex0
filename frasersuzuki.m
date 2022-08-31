function varargout = frasersuzuki(x,p)
   % y = frasersuzuki(x,p)
   %
   % [y, dydp, dydx] = frasersuzuki(x,p)
   %
   % Fraser-Suzuki function (log-normal for r=2, skewed gaussian)
   % See DiMarco2011 - Mathematical functions for the representation of chromatographic peaks.
   %
   % INPUT:      x --> evaluation points
   %             p --> parameter vector where
   %                   p(1) = r  --> narrowness of peak (r=2 for lognormal)
   %                   p(2) = h  --> maximum height of peak
   %                   p(3) = z  --> position of peak
   %                   p(4) = wr --> narrowness of peak
   %                   p(5) = sr --> skewness to the left
   %            
   %
   % OUTPUT:     y --> y(x) values at position x
   %          dydp --> derivative w.r.t. to p
   %                   dy(1) = dydr  --> derivate w.r.t. h
   %                   dy(2) = dydh  --> derivate w.r.t. h
   %                   dy(3) = dydz  --> derivate w.r.t. z
   %                   dy(4) = dydwr --> derivate w.r.t. wr
   %                   dy(5) = dydsr --> derivate w.r.t. sr
   %          dydx --> derivative w.r.t. to x
   %                   dy(6) = dydx  --> derivate w.r.t. x
   %
   % CONSTRAINTS:  sr > 0,  sr != 1,  1 < r < inf, often: r=2.
   %               No check is done on constraints!
   %
   % Andreas Sommer, Aug2017, Aug2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu

   persistent sym_dfdx sym_dfdr sym_dfdh sym_dfdz sym_dfdwr sym_dfdsr sym_f
   persistent fun_dfdx fun_dfdr fun_dfdh fun_dfdz fun_dfdwr fun_dfdsr fun_f
   
   % Generate derivatives ones
   if isempty(sym_dfdx) || (ischar(x) && strcmpi(x,'init'))
      fprintf('frasersuzuki.m: Initializing partial derivatives... ');
      syms sym_r sym_h sym_z sym_wr sym_sr sym_x 
      sym_f = symfun( sym_h .* exp(-(log(sym_r))./(log(sym_sr).^2) * log(((sym_x-sym_z).*(sym_sr.^2-1))./(sym_wr.*sym_sr) + 1).^2)  , ...
                      [sym_x sym_r sym_h sym_z sym_wr sym_sr]) ;  % order of variables
      sym_dfdr  = diff(sym_f,sym_r);
      sym_dfdh  = diff(sym_f,sym_h);
      sym_dfdz  = diff(sym_f,sym_z);
      sym_dfdwr = diff(sym_f,sym_wr);
      sym_dfdsr = diff(sym_f,sym_sr);
      sym_dfdx  = diff(sym_f,sym_x);
      fprintf('Transforming into matlab functions... ');
      fun_f     = matlabFunction(sym_f);
      fun_dfdr  = matlabFunction(sym_dfdr);
      fun_dfdh  = matlabFunction(sym_dfdh);
      fun_dfdz  = matlabFunction(sym_dfdz);
      fun_dfdwr = matlabFunction(sym_dfdwr);
      fun_dfdsr = matlabFunction(sym_dfdsr);
      fun_dfdx  = matlabFunction(sym_dfdx);
      fprintf('Done!\n');
   end

   % accessors
   r  = p(1);
   h  = p(2);
   z  = p(3);
   wr = p(4);
   sr = p(5);
   
   % determine the nonzero indices
   zeropos    = z - (wr.*sr)./(sr^2-1);
   nonzeroIDX = ( x < zeropos );
   xnonzero   = x(nonzeroIDX);
   
   % first output argument: nominal values
   y = zeros(size(x));
   y(nonzeroIDX) = fun_f(xnonzero,r,h,z,wr,sr);
   varargout{1} = y;
   
   % derivatives w.r.t. p requested?
   if (nargout >= 2)
      % derivative w.r.t. r
      dydr = zeros(length(x),1);
      dydr(nonzeroIDX) = fun_dfdr(xnonzero,r,h,z,wr,sr);
      % derivative w.r.t. h
      dydh = zeros(length(x),1);
      dydh(nonzeroIDX) = fun_dfdh(xnonzero,r,h,z,wr,sr);
      % derivative w.r.t. z
      dydz = zeros(length(x),1);
      dydz(nonzeroIDX) = fun_dfdz(xnonzero,r,h,z,wr,sr);
      % derivative w.r.t. wr
      dydwr = zeros(length(x),1);
      dydwr(nonzeroIDX) = fun_dfdwr(xnonzero,r,h,z,wr,sr);
      % derivative w.r.t. sr
      dydsr = zeros(length(x),1);
      dydsr(nonzeroIDX) = fun_dfdsr(xnonzero,r,h,z,wr,sr);
      dydp = [dydr , dydh , dydz , dydwr , dydsr ];
      varargout{2} = dydp;
   end
   % derivatives w.r.t. x requested?
   if (nargout >= 3)
      dydx = zeros(length(x),1);
      dydx(nonzeroIDX) = fun_dfdx(xnonzero,r,h,z,wr,sr);
      varargout{3} = dydx;
   end
  
   
end



% Code for getting the symbolic derivatives
% syms x r h z wr sr
% f = h .* exp(-(log(r))./(log(sr).^2) * log(((x-z)*(sr.^2-1))./(wr.*sr) + 1).^2);
% dfdx  = diff(f,x);
% dfdr  = diff(f,r);
% dfdh  = diff(f,h);
% dfdz  = diff(f,z);
% dfdwr = diff(f,wr);
% dfdsr = diff(f,sr);

