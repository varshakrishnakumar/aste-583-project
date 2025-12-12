function cspice_dasudc( handle, first, last, bpos, epos, data )
   switch nargin
      case 6
         handle = zzmice_int(handle);
         first  = zzmice_int(first);
         last   = zzmice_int(last);
         bpos   = zzmice_int(bpos);
         epos   = zzmice_int(epos);
      otherwise
         error ( [ 'Usage: cspice_dasudc( handle, first, last, bpos, '      ...
                   'epos, data(n,m) )' ] )
   end
   try
      mice('dasudc_c', handle, first, last, bpos, epos, data');
   catch spiceerr
      rethrow(spiceerr)
   end