function [value] = cspice_timdef_get( item )
   switch nargin
      case 1
         item  = zzmice_str(item);
      otherwise
         error ( 'Usage: [value] = cspice_timdef_get( `item` )' )
   end
   try
      [value] = mice( 'timdef_get_c', item );
   catch spiceerr
      rethrow(spiceerr)
   end