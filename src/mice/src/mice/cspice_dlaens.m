function cspice_dlaens( handle )
   switch nargin
      case 1
         handle = zzmice_int(handle);
      otherwise
         error ( 'Usage: cspice_dlaens( handle )' )
   end
   try
      mice('dlaens_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end