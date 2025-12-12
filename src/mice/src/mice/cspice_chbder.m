function [dpdxs] = cspice_chbder( cp, degp, x2s, x, nderiv )
   switch nargin
      case 5
         cp     = zzmice_dp(cp);
         degp   = zzmice_int(degp);
         x2s    = zzmice_dp(x2s);
         x      = zzmice_dp(x);
         nderiv = zzmice_int(nderiv, [0, int32(inf)]);
      otherwise
         error ( [ 'Usage: [dpdxs(nderiv+1)] = '                           ...
                   'cspice_chbder( cp(degp+1), degp, x2s(2), x, nderiv )' ] )
   end
   try
      [dpdxs] = mice('chbder_c', cp, degp, x2s, x, nderiv);
   catch spiceerr
      rethrow(spiceerr)
   end