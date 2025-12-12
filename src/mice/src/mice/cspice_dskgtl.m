function [dpval] = cspice_dskgtl( keywrd )
   switch nargin
      case 1
         keywrd  = zzmice_int(keywrd, [1,6]);
      otherwise
         error ( 'Usage: [dpval] = cspice_dskgtl( keywrd )' )
   end
   try
      [dpval] = mice( 'dskgtl_c', keywrd );
   catch spiceerr
      rethrow(spiceerr)
   end