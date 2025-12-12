function [xform] = cspice_eul2xf(eulang, axisa, axisb, axisc)
   switch nargin
      case 4
         eulang= zzmice_dp(eulang);
         axisa = zzmice_int( axisa);
         axisb = zzmice_int( axisb);
         axisc = zzmice_int( axisc);
      otherwise
         error ( ['Usage: [_xform(6,6)_] = ' ...
                 'cspice_eul2xf(_eulang(6)_, axisa, axisb, axisc)'] )
   end
   try
      [xform] = mice('eul2xf_c', eulang, axisa, axisb, axisc);
   catch spiceerr
      rethrow(spiceerr)
   end