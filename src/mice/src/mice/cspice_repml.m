function [out] = cspice_repml( in, marker, value, rtcase )
   switch nargin
      case 4
         in     = zzmice_str(in, true);
         marker = zzmice_str(marker, true);
         value  = zzmice_int(value);
         rtcase = zzmice_str(rtcase);
      otherwise
         error ( [ 'Usage: [`out`] = '                                      ...
                   'cspice_repml( `in`, `marker`, value, `rtcase` )' ] )
   end
   try
      [out] = mice('repml_c', in, marker, value, rtcase);
   catch spiceerr
      rethrow(spiceerr)
   end