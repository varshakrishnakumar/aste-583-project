function [mout] = cspice_rotmat( m1, angle, iaxis )
   switch nargin
      case 3
         m1    = zzmice_dp(m1);
         angle = zzmice_dp(angle);
         iaxis = zzmice_int(iaxis);
      otherwise
         error ( [ 'Usage: [mout(3,3)] = ' ...
                   'cspice_rotmat( m1(3,3), angle, iaxis )' ] )
   end
   try
      [mout] = mice('rotmat_c', m1, angle, iaxis );
   catch spiceerr
      rethrow(spiceerr)
   end