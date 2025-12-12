function [nresvr, nresvc, ncomr,  ncomc,                                   ...
          free,   lastla, lastrc, lastwd] = cspice_dashfs( handle )
   switch nargin
      case 1
         handle = zzmice_int(handle);
      otherwise
         error ( [ 'Usage: [nresvr, nresvc, ncomr, ncomc, free, '          ...
                   'lastla(3), lastrc(3), lastwd(3)] = '                   ...
                   'cspice_dashfs( handle )' ] )
   end
   try
      [nresvr, nresvc, ncomr,  ncomc,                                      ...
       free,   lastla, lastrc, lastwd] = mice('dashfs_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end