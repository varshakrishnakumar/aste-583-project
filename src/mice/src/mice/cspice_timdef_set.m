function cspice_timdef_set( item, value )
   switch nargin
      case 2
         item  = zzmice_str(item);
         value = zzmice_str(value);
      otherwise
         error ( 'Usage: cspice_timdef_set( `item`, `value`)' )
   end
   try
      mice( 'timdef_set_c', item, value );
   catch spiceerr
      rethrow(spiceerr)
   end