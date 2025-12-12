function cspice_tparch( type )
   switch nargin
      case 1
         type = zzmice_str(type);
      otherwise
         error ( 'Usage: cspice_tparch( `type` )' )
   end
   try
      mice('tparch_c', type);
   catch spiceerr
      rethrow(spiceerr)
   end