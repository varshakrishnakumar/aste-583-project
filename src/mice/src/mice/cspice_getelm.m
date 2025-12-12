function [epoch, elems] = cspice_getelm( frstyr, lines )
   switch nargin
      case 2
         frstyr = zzmice_int(frstyr);
         lines  = zzmice_str(lines);
      otherwise
         error ( ['Usage: [epoch, elems(10) ] = ' ...
                 'cspice_getelm( frstyr, `lines(2)` )'] )
   end
   try
      [epoch, elems] = mice( 'getelm_c', frstyr, lines );
   catch spiceerr
      rethrow(spiceerr)
   end