function [ idata, null, found] = cspice_ekgi( selidx, row, elment )
   switch nargin
      case 3
         selidx = zzmice_int(selidx);
         row    = zzmice_int(row);
         elment = zzmice_int(elment);
      otherwise
         error ( [ 'Usage: [ idata, null, found] = ' ...
                   'cspice_ekgi( selidx, row, elment )' ] )
   end
   try
      [ idata, null, found] = mice('ekgi_c', selidx, row, elment );
      null  = zzmice_logical(null);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end