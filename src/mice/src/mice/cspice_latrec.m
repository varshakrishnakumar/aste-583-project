function [rectan] = cspice_latrec(radius, lon, lat)
   switch nargin
      case 3
         radius = zzmice_dp(radius);
         lon    = zzmice_dp(lon);
         lat    = zzmice_dp(lat);
      otherwise
         error ( ['Usage: [_rectan(3)_] = ' ...
                  'cspice_latrec(_radius_, _lon_, _lat_)'] )
   end
   try
      [rectan] = mice('latrec_c', radius, lon, lat);
   catch spiceerr
      rethrow(spiceerr)
   end