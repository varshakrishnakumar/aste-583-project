function [vnorm] = cspice_vnorm(v1)
   switch nargin
      case 1
         v1 = zzmice_dp(v1);
      otherwise
         error ( 'Usage: [_vnorm_] = cspice_vnorm(_v1(3)_)' )
   end
   try
      [vnorm] = mice('vnorm_c', v1);
   catch spiceerr
      rethrow(spiceerr)
   end