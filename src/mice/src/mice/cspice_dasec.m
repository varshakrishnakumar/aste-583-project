function [buffer, done] = cspice_dasec( handle, bufsiz, buffln )
   switch nargin
      case 3
         handle  = zzmice_int(handle);
         bufsiz  = zzmice_int(bufsiz);
         buffln  = zzmice_int(buffln);
      otherwise
         error ( ['Usage: [buffer, done] = ' ...
                          'cspice_dasec( handle, bufsiz, buffln )'] )
   end
   try
      [buffer, done] = mice( 'dasec_c', handle, bufsiz, buffln );
   catch spiceerr
      rethrow(spiceerr)
      done = zzmice_logical(done);
   end