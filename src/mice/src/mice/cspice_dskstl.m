function cspice_dskstl( keywrd, dpval )
   switch nargin
      case 2
         keywrd  = zzmice_int(keywrd, [1,6]);
         dpval   = zzmice_dp(dpval);
      otherwise
         error ( 'Usage: cspice_dskstl( keywrd, dpval )' )
   end
   try
      mice( 'dskstl_c', keywrd, dpval );
   catch spiceerr
      rethrow(spiceerr)
   end