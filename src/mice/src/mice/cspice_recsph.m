function [r, colat, slon] = cspice_recsph(rectan)
   switch nargin
      case 1
         rectan = zzmice_dp(rectan);
      otherwise
         error ( ['Usage: [_r_, _colat_, _slon_] =' ...
                                ' cspice_recsph( _rectan(3)_ )'] )
   end
   try
      [r, colat, slon] = mice('recsph_c',rectan);
   catch spiceerr
      rethrow(spiceerr)
   end