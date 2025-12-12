function [angle3, angle2, angle1] = cspice_m2eul(r, axis3, axis2, axis1)
   switch nargin
      case 4
         r     = zzmice_dp(r);
         axis3 = zzmice_int(axis3);
         axis2 = zzmice_int(axis2);
         axis1 = zzmice_int(axis1);
      otherwise
         error ( ['Usage: [_angle3_, _angle2_, _angle1_]  = '  ...
                  'cspice_m2eul( _r(3,3)_, axis3, axis2, axis1)']  )
   end
   try
      [angle3, angle2, angle1] = mice('m2eul_c', r, axis3, axis2, axis1);
   catch spiceerr
      rethrow(spiceerr)
   end