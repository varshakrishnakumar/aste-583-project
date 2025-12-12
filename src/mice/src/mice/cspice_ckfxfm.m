function [xform, ref, found] = cspice_ckfxfm( inst, et )
   switch nargin
      case 2
         inst = zzmice_int(inst);
         et   = zzmice_dp(et);
      otherwise
         error ( [ 'Usage: [xform(6,6), ref, found] = '                     ...
                   'cspice_ckfxfm( inst, et )' ] )
   end
   try
      [xform, ref, found] = mice('ckfxfm_c', inst, et);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end