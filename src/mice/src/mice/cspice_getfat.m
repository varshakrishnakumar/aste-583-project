function [arch, kertyp] = cspice_getfat( file )
   switch nargin
      case 1
         file = zzmice_str(file);
      otherwise
         error ( 'Usage: [`arch`, `kertyp`] = cspice_getfat( `file` )' )
   end
   try
      [arch, kertyp] = mice('getfat_c', file);
   catch spiceerr
      rethrow(spiceerr)
   end