function [p] = cspice_vproj( a, b)
   switch nargin
      case 2
         a = zzmice_dp(a);
         b = zzmice_dp(b);
      otherwise
         error ( 'Usage: [_p(3)_] = cspice_vproj(_a(3)_, _b(3)_)' )
   end
   try
      [p] = mice('vproj_c', a, b);
   catch spiceerr
      rethrow(spiceerr)
   end