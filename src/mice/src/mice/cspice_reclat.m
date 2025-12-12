function [radius, lon, lat] = cspice_reclat(rectan)
   switch nargin
      case 1
         rectan = zzmice_dp(rectan);
      otherwise
         error ( ['Usage: [_radius_, _lon_, _lat_] = ' ...
                  'cspice_reclat(_rectan(3)_)'] )
   end
   try
      [radius, lon, lat] = mice('reclat_c',rectan);
   catch spiceerr
      rethrow(spiceerr)
   end