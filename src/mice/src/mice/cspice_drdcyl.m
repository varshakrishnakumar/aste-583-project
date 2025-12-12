function [jacobi] = cspice_drdcyl( r, clon, z )
   switch nargin
      case 3
         r    = zzmice_dp(r);
         clon = zzmice_dp(clon);
         z    = zzmice_dp(z);
      otherwise
         error( 'Usage: [_jacobi(3,3)_] = cspice_drdcyl( _r_, _clon_, _z_)' )
   end
   try
      [jacobi] = mice('drdcyl_c', r, clon, z);
   catch spiceerr
      rethrow(spiceerr)
   end