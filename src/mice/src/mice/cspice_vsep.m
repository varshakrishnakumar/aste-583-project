function [vsep] = cspice_vsep(v1, v2)
   switch nargin
      case 2
         v1 = zzmice_dp(v1);
         v2 = zzmice_dp(v2);
      otherwise
         error ( 'Usage: [_vsep_] = cspice_vsep(_v1(3)_, _v2(3)_)' )
   end
   try
      [vsep] = mice('vsep_c',v1, v2);
   catch spiceerr
      rethrow(spiceerr)
   end