function cspice_lmpool( cvals )
   switch nargin
      case 1
         cvals = zzmice_str( cvals);
      otherwise
         error ( 'Usage: cspice_lmpool( _`cvals`_ )' )
   end
   try
      mice('lmpool_c', cvals );
   catch spiceerr
      rethrow(spiceerr)
   end