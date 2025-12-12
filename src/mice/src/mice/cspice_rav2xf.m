function [xform] = cspice_rav2xf(rot, av)
   switch nargin
      case 2
         rot = zzmice_dp(rot);
         av  = zzmice_dp(av);
      otherwise
         error ( ['Usage: [_xform(6,6)_] = '                               ...
                  'cspice_rav2xf(_rot(3,3)_, _av(3)_)'] )
   end
   try
      [xform] = mice('rav2xf_c', rot, av);
   catch spiceerr
      rethrow(spiceerr)
   end