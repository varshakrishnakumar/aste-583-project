function [out] = cspice_repmc( in, marker, value )
   switch nargin
      case 3
         in     = zzmice_str(in, true);
         marker = zzmice_str(marker, true);
         value  = zzmice_str(value, true);
      otherwise
         error ( 'Usage: [`out`] = cspice_repmc( `in`, `marker`, `value` )' )
   end
   try
      [out] = mice('repmc_c', in, marker, value);
   catch spiceerr
      rethrow(spiceerr)
   end