function [jacobi] = cspice_drdsph( r, colat, slon )
   switch nargin
      case 3
         r     = zzmice_dp(r);
         colat = zzmice_dp(colat);
         slon  = zzmice_dp(slon);
      otherwise
         error( ['Usage: [_jacobi(3,3)_] = ' ...
                         'cspice_drdsph( _r_, _colat_, _slon_ )'] )
   end
   try
      [jacobi] = mice('drdsph_c', r, colat, slon);
   catch spiceerr
      rethrow(spiceerr)
   end