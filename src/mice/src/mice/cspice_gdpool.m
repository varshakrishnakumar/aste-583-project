function [values, found] = cspice_gdpool( name, start, room )
   switch nargin
      case 3
         name  = zzmice_str(name);
         start = zzmice_int(start);
         room  = zzmice_int(room);
      otherwise
         error ( ['Usage: [values(), found] = ' ...
                  'cspice_gdpool( `name`, start, room)' ] )
   end
   try
      [values, found] = mice( 'gdpool_c', name, start, room);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end