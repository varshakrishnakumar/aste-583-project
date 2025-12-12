function [data] = cspice_dasrdc( handle, first, last, bpos, epos, data_i )
   switch nargin
      case 6
         handle = zzmice_int(handle);
         first  = zzmice_int(first);
         last   = zzmice_int(last);
         bpos   = zzmice_int(bpos);
         epos   = zzmice_int(epos);
      otherwise
         error ( [ 'Usage: [data(n,m)] = '                                  ...
                   'cspice_dasrdc( handle, first, last, '                   ...
                   'bpos, epos, data(n,m)' ] )
   end
   try
      [data] = mice('dasrdc_c', handle, first, last, bpos, epos, data_i')';
   catch spiceerr
      rethrow(spiceerr)
   end