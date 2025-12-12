function [idcode] = cspice_ckmeta( ckid, meta )
   switch nargin
      case 2
         ckid = zzmice_int(ckid);
         meta = zzmice_str(meta);
      otherwise
         error ( 'Usage: [idcode] = cspice_ckmeta( ckid, `meta` )' )
   end
   try
      [idcode] = mice('ckmeta_c', ckid, meta);
   catch spiceerr
      rethrow(spiceerr)
   end