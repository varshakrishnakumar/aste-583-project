function cspice_pcpool( name, cvals )
   switch nargin
      case 2
         name  = zzmice_str(name);
         cvals = zzmice_str(cvals);
      otherwise
         error ( 'Usage: cspice_pcpool( `name`, `cvals(n)` )' )
   end
   try
      mice('pcpool_c', name, cvals );
   catch spiceerr
      rethrow(spiceerr)
   end