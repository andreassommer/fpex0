function sigfun = getSigmoidalBaselevelFunction(T,X,varargin)
% sigfun = getSigmoidalBaselevelFunction(T,M,varargin)
% Returns a function delivering the sigmoidal base-level.
%
% INPUT:   T:  x-coordinates (temperatures, times)
%          X:  y-coordinates (1D profile of Cp)
%
% OUTPUT:  sigfun(t): sigmoidal base function
%
% Assumptions:
%   * input consists of (approx.) linear part, curved data, and again a linear part like this:
%                __                                  _
%               /  \                                / \
%              /   _\_______       ===>            /   \
%      _______/___/                          _____/     \___
%
% Algorithm:
%   1) determine slopes of linear parts at begin and end
%   2) generate sigmoidal
%   3) connect the parts
%
%
% Copyright 2016-2022, Andreas Sommer  code@andreas-sommer.eu
%
% Copying and distribution of this file, with or without modification, are permitted in any medium without royalty, 
% provided the copyright notice and this notice are preserved. This file is offered as-is, without any warranty.

   

   lin1start = T(1);   % start of first linear part
   lin1end   = [];     % end   of first linear part
   lin2start = [];     % start of second linear part
   lin2end   = T(end); % end   of second linear part
   midpoint  = [];     % for sigmoidal
   steepness = 1.5;    % for sigmoidal
   
   % process input args
   if hasOption(varargin, 'lin1start'), lin1start = getOption(varargin, 'lin1start'); end
   if hasOption(varargin, 'lin1end'),   lin1end   = getOption(varargin, 'lin1end');   end
   if hasOption(varargin, 'lin2start'), lin2start = getOption(varargin, 'lin2start'); end
   if hasOption(varargin, 'lin2end'),   lin2end   = getOption(varargin, 'lin2end');   end
   if hasOption(varargin, 'steepness'), steepness = getOption(varargin, 'steepness'); end
   if hasOption(varargin, 'midpoint'),  midpoint  = getOption(varargin, 'midpoint');  end

   % assert necessary information is given
   if isempty(lin1end) || isempty(lin2start) 
      error('Must specify lin1end and lin2start !');
   end

   
   % sigmoidal midpoint (if not specified)
   if isempty(midpoint)
      midpoint  = lin1end + (lin2start - lin1end) / 2;
   end

   
   
   % Fit linear parts indices
   lin1idx = (T>=lin1start) & (T<=lin1end);
   lin2idx = (T>=lin2start) & (T<=lin2end);
   lin1p = polyfit(T(lin1idx), X(lin1idx), 1);  % fit linear function
   lin2p = polyfit(T(lin2idx), X(lin2idx), 1);  % to the data
     
   % base level function
   sigfun = @(T) assembleBC(T);
   
   % Finito
   return
   
   
   
   % function for evaluation
   function val = assembleBC(T)
      idx_lin1 = (T<=lin1end);
      idx_lin2 = (T>=lin2start);
      idx_sig  = ~or(idx_lin1, idx_lin2);
      val_lin1 = polyval(lin1p, T(idx_lin1));
      val_lin2 = polyval(lin2p, T(idx_lin2));

      % bleding using sigmoidal
      val_lin1onsig = polyval(lin1p, T(idx_sig));
      val_lin2onsig = polyval(lin2p, T(idx_sig));
      sigmoidal  = 1 ./ (1+exp(-steepness*(T(idx_sig)-midpoint)));
      sigStartVal = val_lin1(end);
      sigEndVal   = val_lin2(1);
      
      % blend depending on direction
      if (sigStartVal < sigEndVal)
         val_sig = val_lin1onsig .* sigmoidal + val_lin2onsig .* (1-sigmoidal);
      else
         val_sig = val_lin1onsig .* (1-sigmoidal) + val_lin2onsig .* sigmoidal;
      end
      
      % assemble and reshape
      val = [val_lin1; val_sig; val_lin2];
      val = reshape(val, size(T));
   end
   
   
end