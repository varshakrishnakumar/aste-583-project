function [jacobi] = cspice_dcyldr( x, y, z)
   switch nargin
      case 3
         x = zzmice_dp(x);
         y = zzmice_dp(y);
         z = zzmice_dp(z);
      otherwise
         error( 'Usage: [_jacobi(3,3)_] = cspice_dcyldr( _x_, _y_, _z_)' )
   end
   try
      [jacobi] = mice('dcyldr_c', x, y, z);
   catch spiceerr
      rethrow(spiceerr)
   end