function [xform] = cspice_sxform(from, to, et)
   switch nargin
      case 3
         from = zzmice_str(from);
         to   = zzmice_str(to);
         et   = zzmice_dp(et);
      otherwise
         error( ['Usage: [_xform(6,6)_] = cspice_sxform( ',           ...
                                              '`from`, `to`, _et_ )'] )
   end
   try
      [xform] = mice('sxform_c', from, to, et);
   catch spiceerr
      rethrow(spiceerr)
   end