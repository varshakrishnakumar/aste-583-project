function [ file, filtyp, srcfil, handle, found ] = cspice_kdata( which, kind )
   switch nargin
      case 2
         which = zzmice_int(which);
         kind  = zzmice_str(kind);
      otherwise
         error( [ 'Usage: [ `file`, `filtyp`, `srcfil`, handle, found ] = ' ...
                                         'cspice_kdata( which, `kind` )']  )
   end
   try
      [ file, filtyp, srcfil, handle, found ]  = mice('kdata_c', which, kind);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end