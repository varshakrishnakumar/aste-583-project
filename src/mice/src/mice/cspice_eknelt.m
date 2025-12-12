function [nelt] = cspice_eknelt( selidx, row )
   switch nargin
      case 2
         selidx = zzmice_int(selidx);
         row    = zzmice_int(row);
      otherwise
         error ( 'Usage: [nelt] = cspice_eknelt( selidx, row )' )
   end
   try
      [nelt] = mice('eknelt_c', selidx, row );
   catch spiceerr
      rethrow(spiceerr)
   end