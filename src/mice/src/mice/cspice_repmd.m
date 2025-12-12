function [out] = cspice_repmd( in, marker, value, sigdig )
   switch nargin
      case 4
         in     = zzmice_str(in, true);
         marker = zzmice_str(marker, true);
         value  = zzmice_dp(value);
         sigdig = zzmice_int(sigdig);
      otherwise
         error ( [ 'Usage: [`out`] = '                                      ...
                   'cspice_repmd( `in`, `marker`, value, sigdig )' ] )
   end
   try
      [out] = mice('repmd_c', in, marker, value, sigdig);
   catch spiceerr
      rethrow(spiceerr)
   end