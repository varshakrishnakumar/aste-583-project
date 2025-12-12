function [p] = cspice_chbval( cp, degp, x2s, x )
   switch nargin
      case 4
         cp   = zzmice_dp(cp);
         degp = zzmice_int(degp);
         x2s  = zzmice_dp(x2s);
         x    = zzmice_dp(x);
      otherwise
         error ( 'Usage: [p] = cspice_chbval( cp(degp+1), degp, x2s(2), x )' )
   end
   try
      [p] = mice('chbval_c', cp, degp, x2s, x);
   catch spiceerr
      rethrow(spiceerr)
   end