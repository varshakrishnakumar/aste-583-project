function [radius, colat, slon] = cspice_cylsph( r, clon, z )
   switch nargin
      case 3
         r    = zzmice_dp(r);
         clon = zzmice_dp(clon);
         z    = zzmice_dp(z);
      otherwise
         error ( ['Usage: [_radius_, _colat_, _slon_] = '...
                                     'cspice_cylsph(_r_, _clon_, _z_)' ] )
   end
   try
      [radius, colat, slon] = mice('cylsph_c', r, clon, z);
   catch spiceerr
      rethrow(spiceerr)
   end