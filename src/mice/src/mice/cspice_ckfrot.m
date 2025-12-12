function [rotate, ref, found] = cspice_ckfrot( inst, et )
   switch nargin
      case 2
         inst = zzmice_int(inst);
         et   = zzmice_dp(et);
      otherwise
         error ( [ 'Usage: [rotate(3,3), ref, found] = '                    ...
                   'cspice_ckfrot( inst, et )' ] )
   end
   try
      [rotate, ref, found] = mice('ckfrot_c', inst, et);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end