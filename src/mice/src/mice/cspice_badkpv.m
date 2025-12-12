function cspice_badkpv( caller, name, comp, size, divby, type )
   switch nargin
      case 6
         caller = zzmice_str(caller);
         name   = zzmice_str(name);
         comp   = zzmice_str(comp);
         size   = zzmice_int(size);
         divby  = zzmice_int(divby);
         type   = zzmice_str(type);
      otherwise
         error ( [ 'Usage: '                                                ...
                   'cspice_badkpv( `caller`, `name`, `comp`, size, divby, ' ...
                   '`type` )' ] )
   end
   try
      mice('badkpv_c', caller, name, comp, size, divby, type);
   catch spiceerr
      rethrow(spiceerr)
   end