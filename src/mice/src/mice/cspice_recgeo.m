function [lon, lat, alt] = cspice_recgeo(rectan, re, f)
   switch nargin
      case 3
         rectan = zzmice_dp(rectan);
         re     = zzmice_dp(re);
         f      = zzmice_dp(f);
      otherwise
         error ( ['Usage: [_lon_, _lat_, _alt_] = '...
                  'cspice_recgeo(_rectan(3)_, re, f)' ] )
   end
   try
      [lon, lat, alt] = mice( 'recgeo_c', rectan, re, f);
   catch spiceerr
      rethrow(spiceerr)
   end