function [eulang, unique] = cspice_xf2eul(xform,  axisa, axisb, axisc)
   switch nargin
      case 4
         xform = zzmice_dp( xform);
         axisa = zzmice_int( axisa);
         axisb = zzmice_int( axisb);
         axisc = zzmice_int( axisc);
      otherwise
         error ( [ 'Usage: [_eulang(6)_, _unique_] = ' ...
                   'cspice_xf2eul(_xform(6,6)_, axisa, axisb, axisc)'] )
   end
   try
      [eulang, unique] = mice('xf2eul_c', xform, axisa, axisb, axisc);
      unique = zzmice_logical(unique);
   catch spiceerr
      rethrow(spiceerr)
   end