function [found, n, type] = cspice_dtpool(name)
   switch nargin
      case 1
         name = zzmice_str(name);
      otherwise
         error ( 'Usage: [_found_, _n_, _`type`_] = cspice_dtpool(_`name`_)' )
   end
   try
      [dtpool] = mice('dtpool_s',name);
      found    = reshape( [dtpool.found], 1, [] );
      n        = reshape( [dtpool.n],     1, [] );
      type     = char( dtpool.type );
   catch spiceerr
      rethrow(spiceerr)
   end