function [radius, lon, lat] = cspice_cyllat(r, clon, z)
   switch nargin
      case 3
         r    = zzmice_dp(r);
         clon = zzmice_dp(clon);
         z    = zzmice_dp(z);
      otherwise
         error ( ['Usage: [_radius_, _lon_, _lat_] = '...
                  'cspice_cyllat(_r_, _clon_, _z_)'] )
   end
   try
      [radius, lon, lat] = mice('cyllat_c', r, clon, z);
   catch spiceerr
      rethrow(spiceerr)
   end