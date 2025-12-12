function [verstr] = cspice_tkvrsn( item )
   switch nargin
      case 1
         item = zzmice_str(item);
      otherwise
         error ( 'Usage: [`verstr`] = cspice_tkvrsn( `item` )' )
   end
   try
      [verstr] = mice( 'tkvrsn_c', item );
   catch spiceerr
      rethrow(spiceerr)
   end