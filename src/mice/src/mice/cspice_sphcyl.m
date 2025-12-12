function [ r, clon, z] = cspice_sphcyl(radius, colat, slon)
   switch nargin
      case 3
         radius = zzmice_dp(radius);
         colat  = zzmice_dp(colat);
         slon   = zzmice_dp(slon);
      otherwise
         error ( ['Usage: [ _r_, _clon_, _z_] = '...
                  'cspice_sphcyl(_radius_, _colat_, _slon_)' ] )
   end
   try
      [ r, clon, z] = mice('sphcyl_c', radius, colat, slon );
   catch spiceerr
      rethrow(spiceerr)
   end