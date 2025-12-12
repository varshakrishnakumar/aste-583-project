function [ep] = cspice_edpnt( p, a, b, c )
   switch nargin
      case 4
         p = zzmice_dp(p);
         a = zzmice_dp(a);
         b = zzmice_dp(b);
         c = zzmice_dp(c);
      otherwise
         error ( 'Usage: [ep(3)] = cspice_edpnt( p(3), a, b, c )' )
   end
   try
      [ep] = mice('edpnt_c', p, a, b, c);
   catch spiceerr
      rethrow(spiceerr)
   end