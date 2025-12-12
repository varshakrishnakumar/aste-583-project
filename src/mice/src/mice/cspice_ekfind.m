function [nmrows, error, errmsg] = cspice_ekfind(query)
   switch nargin
      case 1
         sample = zzmice_str(query);
      otherwise
         error ( [ 'Usage: [nmrows, error, `errmsg`] = ' ...
                   'cspice_ekfind( `query` )' ])
   end
   try
      [nmrows, error, errmsg] = mice('ekfind_c', query );
      error = zzmice_logical(error);
   catch spiceerr
      rethrow(spiceerr)
   end