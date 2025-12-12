function [vout] = cspice_vhat(v1)
   switch nargin
      case 1
         v1 = zzmice_dp(v1);
      otherwise
         error ( 'Usage: [_vout(3)_] = cspice_vhat(_v1(3)_)' )
   end
   try
      [vout] = mice('vhat_c',v1);
   catch spiceerr
      rethrow(spiceerr)
   end