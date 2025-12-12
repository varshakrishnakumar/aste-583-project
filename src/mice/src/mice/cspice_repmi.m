function [out] = cspice_repmi( in, marker, value )
   switch nargin
      case 3
         in     = zzmice_str(in, true);
         marker = zzmice_str(marker, true);
         value  = zzmice_int(value);
      otherwise
         error ( 'Usage: [`out`] = cspice_repmi( `in`, `marker`, value )' )
   end
   try
      [out] = mice('repmi_c', in, marker, value);
   catch spiceerr
      rethrow(spiceerr)
   end