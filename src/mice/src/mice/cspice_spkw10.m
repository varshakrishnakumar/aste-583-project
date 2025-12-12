function cspice_spkw10( handle, body,   center, frame, first,  ...
                        last,   segid,  consts, n,     elems,  ...
                        epochs )
   switch nargin
      case 11
         handle = zzmice_int(handle);
         body   = zzmice_int(body);
         center = zzmice_int(center);
         frame  = zzmice_str(frame);
         first  = zzmice_dp(first);
         last   = zzmice_dp(last);
         segid  = zzmice_str(segid);
         consts = zzmice_dp(consts);
         n      = zzmice_int(n);
         elems  = zzmice_dp(elems);
         epochs = zzmice_dp(epochs);
      otherwise
         error ( [ 'Usage: cspice_spkw10( handle, body, center, `frame`, '  ...
                   'first, last, `segid`, consts(8), n, elems(10*n), '      ...
                   'epochs(n) )=' ] )
   end
   try
      mice( 'spkw10_c', handle, body,   center, frame,  first,              ...
            last,       segid,  consts, n,      elems,  epochs );
   catch spiceerr
      rethrow(spiceerr)
   end