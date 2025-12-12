function [vout] = cspice_rotvec( v1, angle, iaxis )
   switch nargin
      case 3
         v1 = zzmice_dp(v1);
         angle = zzmice_dp(angle);
         iaxis = zzmice_int(iaxis);
      otherwise
         error ( 'Usage: [vout(3)] = cspice_rotvec( v1(3), angle, iaxis )' )
   end
   try
      [vout] = mice('rotvec_c', v1, angle, iaxis);
   catch spiceerr
      rethrow(spiceerr)
   end