function [jacobi] = cspice_drdgeo( lon, lat, alt, re, f )
   switch nargin
      case 5
         lon = zzmice_dp(lon);
         lat = zzmice_dp(lat);
         alt = zzmice_dp(alt);
         re  = zzmice_dp(re);
         f   = zzmice_dp(f);
      otherwise
         error( ['Usage: [_jacobi(3,3)_] = '...
                 'cspice_drdgeo( _lon_, _lat_, _alt_, re, f)'] )
   end
   try
      [jacobi] = mice('drdgeo_c', lon, lat, alt, re, f);
   catch spiceerr
      rethrow(spiceerr)
   end