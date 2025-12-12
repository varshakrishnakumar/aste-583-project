function [lon, lat, alt] = cspice_recpgr( body, rectan, re, f)
   switch nargin
      case 4
         body   = zzmice_str(body);
         rectan = zzmice_dp(rectan);
         re     = zzmice_dp(re);
         f      = zzmice_dp(f);
      otherwise
         error ( ['Usage: [_lon_, _lat_, _alt_] = '...
                  'cspice_recpgr( body, _rectan(3)_, re, f)' ] )
   end
   try
      [lon, lat, alt] = mice( 'recpgr_c', body, rectan, re, f);
   catch spiceerr
      rethrow(spiceerr)
   end