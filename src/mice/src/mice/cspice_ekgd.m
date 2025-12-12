function [ ddata, null, found] = cspice_ekgd( selidx, row, elment )
   switch nargin
      case 3
         selidx = zzmice_int(selidx);
         row    = zzmice_int(row);
         elment = zzmice_int(elment);
      otherwise
         error ( [ 'Usage: [ ddata, null, found] = ' ...
                   'cspice_ekgd( selidx, row, elment )' ] )
   end
   try
      [ ddata, null, found] = mice('ekgd_c', selidx, row, elment );
      null  = zzmice_logical(null);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end