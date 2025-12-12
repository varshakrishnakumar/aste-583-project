function [rectan] = cspice_cylrec(r, clon, z)
   switch nargin
      case 3
         r    = zzmice_dp(r);
         clon = zzmice_dp(clon);
         z    = zzmice_dp(z);
      otherwise
         error ( 'Usage: [_rectan(3)_] = cspice_cylrec(_r_, _clon_, _z_)' )
   end
   try
      [rectan] = mice('cylrec_c', r, clon, z);
   catch spiceerr
      rethrow(spiceerr)
   end