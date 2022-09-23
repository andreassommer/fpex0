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
   % This is a symbolic-free implementation. See the frasersuzuki_symbolic.m file to learn
   % how to easily generate initial distribution functions with automatic derivative generation
   % if Matlab's Symbolic Math Toolbox is available.
   %
   % Andreas Sommer, Aug2017, Aug2022, Sep2022
   % andreas.sommer@iwr.uni-heidelberg.de
   % code@andreas-sommer.eu


 
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
   

   
   % common subexpressions (for f and df)
   Xmz        = (xnonzero-z);
   sr2m1      = (sr^2-1);
   Xmz_sr2m1_by_wrsr_p1 = (Xmz*sr2m1)/(wr*sr) + 1;
   log_X      = log( Xmz_sr2m1_by_wrsr_p1 );
   log_X_sq   = log_X .^ 2;
   log_sr2    = log(sr)^2;
   log_r_by_log_sr2 = log(r) / log_sr2;
   exp_X      = exp(- log_X_sq * log_r_by_log_sr2);
   h_exp_X    = h * exp_X;


   % FIRST OUTPUT: nominal values
   y = zeros(size(x));
   y(nonzeroIDX) = h_exp_X;
   varargout{1} = y;


   % derivatives w.r.t. p requested?
   if (nargout >= 2)

      % common subexpressions (for df only)
      h_exp_X_log_X = h_exp_X .* log_X;
      h_exp_X_log_X_sq = h_exp_X .* log_X_sq;
      log_sr2_X = log(sr)^2 * Xmz_sr2m1_by_wrsr_p1;
      
      % derivates w.r.t. parameters (computed only at non-zero values)
      dfdr  = - (h_exp_X_log_X_sq ) / (r*log(sr)^2);
      dfdh  = exp_X;
      dfdz  = (2*h_exp_X_log_X * log(r)*sr2m1)          ./ (sr*wr * log_sr2_X );
      dfdwr = (2*h_exp_X_log_X * log(r)*sr2m1 .* Xmz)   ./ (sr*wr^2*log_sr2_X );
      dfdsr = h_exp_X .* (   (2*log_X_sq * log(r)) / (sr*log(sr)^3)    ...
                           - (2*log_X    * log(r) .* ( 2*Xmz/wr - (sr2m1 * Xmz) / (sr^2*wr))) ./ log_sr2_X  );
   
      % transfer into output variables
      dydr  = zeros(length(x),1);  dydr(nonzeroIDX)  = dfdr; 
      dydh  = zeros(length(x),1);  dydh(nonzeroIDX)  = dfdh;
      dydz  = zeros(length(x),1);  dydz(nonzeroIDX)  = dfdz; 
      dydwr = zeros(length(x),1);  dydwr(nonzeroIDX) = dfdwr;
      dydsr = zeros(length(x),1);  dydsr(nonzeroIDX) = dfdsr;
      dydp = [dydr , dydh , dydz , dydwr , dydsr ];
      varargout{2} = dydp;
   end
   
   
   % derivatives w.r.t. x requested?
   if (nargout >= 3)
      % derivates w.r.t. x (computed only at non-zero values)
      dfdx  = -(2*h_exp_X_log_X * log(r)*(sr^2 - 1))    ./ (sr*wr * log_sr2_X );
      % transfer into output variables
      dydx = zeros(length(x),1);  dydx(nonzeroIDX) = dfdx;
      varargout{3} = dydx;
   end

   
   
   % finito
   return
   
   
   
   %% Helper functions (nominal value and derivatives, not used)
%    function val = fun_f(x)
%       val =  h * exp(-(log(r))/(log(sr)^2) * log( ((x-z)*(sr^2-1))/(wr*sr) + 1).^2 );
%    end
%    function val = fun_dfdr(x)
%       val = -(h*exp(-(log(((sr^2-1)*(x-z))/(sr*wr) + 1).^2*log(r))/log(sr)^2)  ...
%                .*log(((sr^2-1)*(x-z))/(sr*wr)+1).^2) / (r*log(sr)^2);
%    end
%    function val = fun_dfdh(x)
%       val = exp( -(log(((sr^2 - 1)*(x - z))/(sr*wr) + 1).^2*log(r)) / log(sr)^2 );
%    end
%    function val = fun_dfdz(x)
%       val = (2*h*exp(-(log(((sr^2 - 1)*(x - z))/(sr*wr) + 1).^2*log(r))/log(sr)^2)   ...
%                .*log(((sr^2 - 1)*(x - z))/(sr*wr) + 1)                              ...
%                 *log(r)*(sr^2 - 1))                                                 ...
%             ./ (sr*wr*log(sr)^2 * (((sr^2 - 1)*(x - z))/(sr*wr) + 1));
%    end
%    function val = fun_dfdwr(x)
%       val = (2*h*exp(-(log(((sr^2 - 1)*(x - z))/(sr*wr) + 1).^2*log(r)) / log(sr)^2) ...
%                .*log(((sr^2 - 1)*(x - z))/(sr*wr) + 1)                              ...
%                 *log(r)*(sr^2 - 1)                                                  ...
%                .*(x - z))                                                           ...
%             ./(sr*wr^2*log(sr)^2*(((sr^2 - 1)*(x - z))/(sr*wr) + 1));
%    end
%    function val = fun_dfdsr(x)
%       val = h*exp(-(log(((sr^2 - 1)*(x - z))/(sr*wr) + 1).^2*log(r))/log(sr)^2)               ...
%             .*( (2*log(((sr^2 - 1)*(x - z))/(sr*wr) + 1).^2 * log(r)) / (sr*log(sr)^3)       ...
%                 - (2*log(((sr^2 - 1)*(x - z))/(sr*wr) + 1)                       ...
%                    .*log(r).*((2*(x - z))/wr - ((sr^2 - 1)*(x - z))/(sr^2*wr)))   ...
%                    ./ ( log(sr)^2*(((sr^2 - 1)*(x - z))/(sr*wr) + 1))            ...
%               );
%    end
%    function val = fun_dfdx(x)
%       val = -(2*h*exp(-(log(((sr^2 - 1)*(x - z))/(sr*wr) + 1).^2*log(r))/log(sr)^2)*log(((sr^2 - 1)*(x - z))/(sr*wr) + 1)*log(r)*(sr^2 - 1))   ...
%              ./ (sr*wr*log(sr)^2*(((sr^2 - 1)*(x - z))/(sr*wr) + 1));
%    end
      
   
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

