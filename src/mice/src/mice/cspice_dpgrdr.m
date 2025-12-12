function [jacobi] = cspice_dpgrdr( body, x, y, z, re, f)
   switch nargin
      case 6
         body = zzmice_str(body);
         x    = zzmice_dp(x);
         y    = zzmice_dp(y);
         z    = zzmice_dp(z);
         re   = zzmice_dp(re);
         f    = zzmice_dp(f);
      otherwise
         error( ['Usage: [_jacobi(3,3)_] = ' ...
                 'cspice_dpgrdr( `body`, _x_, _y_, _z_, re, f)'] )
   end
   try
      [jacobi] = mice('dpgrdr_c', body, x, y, z, re, f);
   catch spiceerr
      rethrow(spiceerr)
   end