function [p] = cspice_vprojg( a, b )
   switch nargin
      case 2
         a = zzmice_dp(a);
         b = zzmice_dp(b);
      otherwise
         error ( 'Usage: [p(n)] = cspice_vprojg( a(n), b(n) )' )
   end
   try
      [p] = mice('vprojg_c', a, b);
   catch spiceerr
      rethrow(spiceerr)
   end