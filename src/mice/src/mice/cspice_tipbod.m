function [tipm] = cspice_tipbod( ref, body, et )
   switch nargin
      case 3
         ref  = zzmice_str(ref);
         body = zzmice_int(body);
         et   = zzmice_dp(et);
      otherwise
         error ( 'Usage: [tipm(3,3)] = cspice_tipbod( `ref`, body, et )' )
   end
   try
      [tipm] = mice('tipbod_c', ref, body, et);
   catch spiceerr
      rethrow(spiceerr)
   end