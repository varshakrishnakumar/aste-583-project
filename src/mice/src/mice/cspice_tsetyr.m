function cspice_tsetyr( year )
   switch nargin
      case 1
         year = zzmice_int(year);
      otherwise
         error ( 'Usage: cspice_tsetyr( year )' )
   end
   try
      mice('tsetyr_c', year);
   catch spiceerr
      rethrow(spiceerr)
   end