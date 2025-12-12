function [cdata, null, found] = cspice_ekgc( selidx, row, elment, cdatln )
   switch nargin
      case 4
         selidx = zzmice_int(selidx);
         row    = zzmice_int(row);
         elment = zzmice_int(elment);
         cdatln = zzmice_int(cdatln);
      otherwise
         error ( [ 'Usage: [ `cdata`, null, found] = ' ...
                   'cspice_ekgc( selidx, row, elment, cdatln )' ] )
   end
   try
      [cdata, null, found] = mice('ekgc_c', selidx, row, elment, cdatln );
      null  = zzmice_logical(null);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end