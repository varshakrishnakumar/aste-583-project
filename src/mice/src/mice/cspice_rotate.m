function [mout] = cspice_rotate( angle, iaxis )
   switch nargin
      case 2
         angle = zzmice_dp(angle);
         iaxis = zzmice_int( iaxis );
      otherwise
         error ( 'Usage: [_mout(3,3)_] = cspice_rotate( _angle_, iaxis )' )
   end
   try
      [mout] = mice('rotate_c', angle, iaxis );
   catch spiceerr
      rethrow(spiceerr)
   end