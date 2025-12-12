function [radius, lon, lat] = cspice_sphlat(r, colat, slon)
   switch nargin
      case 3
         r     = zzmice_dp(r);
         colat = zzmice_dp(colat);
         slon  = zzmice_dp(slon);
      otherwise
         error ( [ 'Usage: [_radius_, _lon_, _lat_] = ' ...
                   'cspice_sphlat(_r_, _colat_, _slon_)' ] )
   end
   try
      [radius, lon, lat] = mice('sphlat_c', r, colat, slon);
   catch spiceerr
      rethrow(spiceerr)
   end