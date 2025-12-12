function [jacobi] = cspice_dgeodr( x, y, z, re, f)
   switch nargin
      case 5
         x = zzmice_dp(x);
         y = zzmice_dp(y);
         z = zzmice_dp(z);
         re= zzmice_dp(re);
         f = zzmice_dp(f);
      otherwise
         error( [ 'Usage: [_jacobi(3,3)_] = '                              ...
                  'cspice_dgeodr( _x_, _y_, _z_, re, f )' ])
   end
   try
      [jacobi] = mice('dgeodr_c', x, y, z, re, f);
   catch spiceerr
      rethrow(spiceerr)
   end