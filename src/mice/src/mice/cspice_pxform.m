function [rotate] = cspice_pxform(from, to, et)
   switch nargin
      case 3
         from = zzmice_str(from);
         to   = zzmice_str(to);
         et   = zzmice_dp(et);
      otherwise
         error ( [ 'Usage: [_rotate(3,3)_] = ' ...
                   'cspice_pxform( `from`, `to`, _et_ )' ] )
   end
   try
      [rotate] = mice('pxform_c',from,to,et);
   catch spiceerr
      rethrow(spiceerr)
   end