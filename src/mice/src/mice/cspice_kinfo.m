function [filtyp, srcfil, handle, found] = cspice_kinfo( file )
   switch nargin
      case 1
         file = zzmice_str(file);
      otherwise
         error( [ 'Usage: [ `filtyp`, `srcfil`, handle, found ] = ' ...
                                             'cspice_kinfo( `file` )']  )
   end
   try
      [filtyp, srcfil, handle, found]  = mice('kinfo_c', file);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end