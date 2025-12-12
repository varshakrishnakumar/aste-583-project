function [rectan] = cspice_sphrec( r, colat, slon )
   switch nargin
      case 3
         r     = zzmice_dp(r);
         colat = zzmice_dp(colat);
         slon  = zzmice_dp(slon);
      otherwise
         error ( ['Usage: [_rectan(3)_] = ' ...
                        'cspice_sphrec( _r_, _colat_, _slon_ )'] )
   end
   try
      [rectan] = mice('sphrec_c', r, colat, slon);
   catch spiceerr
      rethrow(spiceerr)
   end