function [r] = cspice_eul2m(angle3, angle2, angle1, axis3, axis2, axis1)
   switch nargin
      case 6
         angle3 = zzmice_dp(angle3);
         angle2 = zzmice_dp(angle2);
         angle1 = zzmice_dp(angle1);
         axis3  = zzmice_int(axis3);
         axis2  = zzmice_int(axis2);
         axis1  = zzmice_int(axis1);
      otherwise
         error( ['Usage: [_r(3,3)_] = '                         ...
                 'cspice_eul2m(_angle3_, _angle2_, _angle1_, '  ...
                 'axis3, axis2, axis1)']  )
   end
   try
      [r] = mice('eul2m_c',angle3,angle2,angle1,axis3,axis2,axis1);
   catch spiceerr
      rethrow(spiceerr)
   end