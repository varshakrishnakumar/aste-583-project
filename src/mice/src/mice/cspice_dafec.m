function [buffer, done] = cspice_dafec( handle, bufsiz, buffln )
   switch nargin
      case 3
         handle  = zzmice_int(handle);
         bufsiz  = zzmice_int(bufsiz);
         buffln  = zzmice_int(buffln);
      otherwise
         error ( ['Usage: [buffer, done] = ' ...
                          'cspice_dafec( handle, bufsiz, buffln )'] )
   end
   try
      [buffer, done] = mice( 'dafec_c', handle, bufsiz, buffln );
      done = zzmice_logical(done);
   catch spiceerr
      rethrow(spiceerr)
   end