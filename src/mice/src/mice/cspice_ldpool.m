function cspice_ldpool( fname )
   switch nargin
      case 1
         fname = zzmice_str(fname);
      otherwise
         error ( 'Usage: cspice_ldpool( `fname` )' )
   end
   try
      mice('ldpool_c', fname);
   catch spiceerr
      rethrow(spiceerr)
   end