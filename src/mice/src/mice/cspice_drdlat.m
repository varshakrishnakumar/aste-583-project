function [jacobi] = cspice_drdlat( r, lon, lat)
   switch nargin
      case 3
         r   = zzmice_dp(r);
         lon = zzmice_dp(lon);
         lat = zzmice_dp(lat);
      otherwise
         error( 'Usage: [_jacobi(3,3)_] = cspice_drdlat( _r_, _lon_, _lat_)' )
   end
   try
      [jacobi] = mice('drdlat_c', r, lon, lat);
   catch spiceerr
      rethrow(spiceerr)
   end