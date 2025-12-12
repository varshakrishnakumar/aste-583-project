function [jacobi] = cspice_dsphdr( x, y, z)
   switch nargin
      case 3
         x = zzmice_dp(x);
         y = zzmice_dp(y);
         z = zzmice_dp(z);
      otherwise
         error( 'Usage: [_jacobi(3,3)_] = cspice_dsphdr( _x_, _y_, _z_)' )
   end
   try
      [jacobi] = mice('dsphdr_c', x, y, z);
   catch spiceerr
      rethrow(spiceerr)
   end