function [p] = cspice_polyds( coeffs, deg, nderiv, t )
   switch nargin
      case 4
         coeffs = zzmice_dp(coeffs);
         deg    = zzmice_int(deg);
         nderiv = zzmice_int(nderiv, [0, int32(inf)]);
         t      = zzmice_dp(t);
      otherwise
         error ( [ 'Usage: [p(nderiv+1)] = '                               ...
                   'cspice_polyds( coeffs(deg+1), deg, nderiv, t )' ] )
   end
   try
      [p] = mice('polyds_c', coeffs, deg, nderiv, t);
   catch spiceerr
      rethrow(spiceerr)
   end