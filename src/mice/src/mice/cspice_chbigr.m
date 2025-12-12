function [p, itgrlp] = cspice_chbigr( degp, cp, x2s, x )
   switch nargin
      case 4
         degp = zzmice_int(degp);
         cp   = zzmice_dp(cp);
         x2s  = zzmice_dp(x2s);
         x    = zzmice_dp(x);
      otherwise
         error ( [ 'Usage: [p, itgrlp] = '                                 ...
                   'cspice_chbigr( degp, cp(degp+1), x2s(2), x )' ] )
   end
   try
      [p, itgrlp] = mice('chbigr_c', degp, cp, x2s, x);
   catch spiceerr
      rethrow(spiceerr)
   end