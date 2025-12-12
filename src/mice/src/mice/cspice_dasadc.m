function cspice_dasadc( handle, n, bpos, epos, data )
   switch nargin
      case 5
         handle = zzmice_int(handle);
         n      = zzmice_int(n);
         bpos   = zzmice_int(bpos);
         epos   = zzmice_int(epos);
      otherwise
         error ( 'Usage: cspice_dasadc( handle, n, bpos, epos, data(n,m) )' )
   end
   try
      mice('dasadc_c', handle, n, bpos, epos, data');
   catch spiceerr
      rethrow(spiceerr)
   end