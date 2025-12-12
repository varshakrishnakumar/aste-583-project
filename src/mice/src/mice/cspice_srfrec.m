function [rectan] = cspice_srfrec(body, lon, lat)
   switch nargin
      case 3
         body = zzmice_int(body);
         lon  = zzmice_dp(lon);
         lat  = zzmice_dp(lat);
      otherwise
         error ( ['Usage: [_rectan(3)_] = ' ...
                  'cspice_srfrec(body, _lon_, _lat_)'] )
   end
   try
      [rectan] = mice('srfrec_c', body, lon, lat);
   catch spiceerr
      rethrow(spiceerr)
   end