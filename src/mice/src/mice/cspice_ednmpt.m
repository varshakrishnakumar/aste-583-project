function [point] = cspice_ednmpt( a, b, c, normal )
   switch nargin
      case 4
         a = zzmice_dp(a);
         b = zzmice_dp(b);
         c = zzmice_dp(c);
         normal = zzmice_dp(normal);
      otherwise
         error ( 'Usage: [point(3)] = cspice_ednmpt( a, b, c, normal(3) )' )
   end
   try
      [point] = mice('ednmpt_c', a, b, c, normal);
   catch spiceerr
      rethrow(spiceerr)
   end