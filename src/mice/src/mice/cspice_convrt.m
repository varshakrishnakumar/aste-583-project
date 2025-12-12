function [y] = cspice_convrt( x, in, out)
   switch nargin
      case 3
         x   = zzmice_dp(x);
         in  = zzmice_str(in);
         out = zzmice_str(out);
      otherwise
         error( 'Usage: [_y_] = cspice_convrt( _x_, `in`, `out`)' )
   end
   try
      [y] = mice('convrt_c', x, in, out);
   catch spiceerr
      rethrow(spiceerr)
   end