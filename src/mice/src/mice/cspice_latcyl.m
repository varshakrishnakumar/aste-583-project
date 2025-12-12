function [r, clon, z ] = cspice_latcyl( radius, lon, lat)
   switch nargin
      case 3
         radius = zzmice_dp(radius);
         lon    = zzmice_dp(lon);
         lat    = zzmice_dp(lat);
      otherwise
         error ( [ 'Usage: [ _r_, _clon_, _z_] = '...
                   'cspice_latcyl( _radius_, _lon_, _lat_)' ] )
   end
   try
      [r, clon, z] = mice('latcyl_c', radius, lon, lat);
   catch spiceerr
      rethrow(spiceerr)
   end