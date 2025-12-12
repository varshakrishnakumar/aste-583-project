function [nthstr, found] = cspice_stpool( item, nth, contin )
   switch nargin
      case 3
         item   = zzmice_str(item);
         nth    = zzmice_int(nth);
         contin = zzmice_str(contin);
      otherwise
         error ( ['Usage: [`nthstr`, found] = '                            ...
                   'cspice_stpool( `item`, nth, `contin` )' ] )
   end
   try
      [nthstr, found] = mice( 'stpool_c', item, nth, contin );
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end