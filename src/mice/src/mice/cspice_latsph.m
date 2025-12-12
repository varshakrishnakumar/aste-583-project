function [rho, colat, slon] = cspice_latsph( radius, lon, lat)
   switch nargin
      case 3
         radius = zzmice_dp(radius);
         lon    = zzmice_dp(lon);
         lat    = zzmice_dp(lat);
      otherwise
         error ( ['Usage: [_rho_, _colat_, _slon_] = '...
                  'cspice_latsph( _radius_, _lon_, _lat_)' ] )
   end
   try
      [rho, colat, slon] = mice('latsph_c', radius, lon, lat);
   catch spiceerr
      rethrow(spiceerr)
   end