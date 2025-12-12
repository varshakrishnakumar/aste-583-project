function [rot, av] = cspice_xf2rav(xform)
   switch nargin
      case 1
         xform = zzmice_dp(xform);
      otherwise
         error ( ['Usage: [_rot(3,3)_, _av(3)_] = ' ...
                  'cspice_xf2rav(_xform(6,6)_)'] )
   end
   try
      [rot, av] = mice('xf2rav_c', xform);
   catch spiceerr
      rethrow(spiceerr)
   end