function [vout, vmag] = cspice_unorm( v1 )
   switch nargin
      case 1
         v1 = zzmice_dp(v1);
      otherwise
         error ( 'Usage: [_vout(3)_, _vmag_] = cspice_unorm(_v1(3)_)' )
   end
   try
      [vout, vmag] = mice('unorm_c',v1);
   catch spiceerr
      rethrow(spiceerr)
   end