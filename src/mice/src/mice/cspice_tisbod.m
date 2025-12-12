function [tsipm] = cspice_tisbod( ref, body, et )
   switch nargin
      case 3
         ref  = zzmice_str(ref);
         body = zzmice_int(body);
         et   = zzmice_dp(et);
      otherwise
         error ( 'Usage: [tsipm(6,6)] = cspice_tisbod( `ref`, body, et )' )
   end
   try
      [tsipm] = mice('tisbod_c', ref, body, et);
   catch spiceerr
      rethrow(spiceerr)
   end