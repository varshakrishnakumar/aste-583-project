function [rectan] = cspice_pgrrec( body, lon, lat, alt, re, f)
   switch nargin
      case 6
         body = zzmice_str(body);
         lon  = zzmice_dp(lon);
         lat  = zzmice_dp(lat);
         alt  = zzmice_dp(alt);
         re   = zzmice_dp(re);
         f    = zzmice_dp(f);
      otherwise
         error ( ['Usage: [_rectan(3)_] = ' ...
                  'cspice_pgrrec( body, _lon_, _lat_, _alt_, re, f)'] )
   end
   try
      [rectan] = mice( 'pgrrec_c', body, lon, lat, alt, re, f);
   catch spiceerr
      rethrow(spiceerr)
   end