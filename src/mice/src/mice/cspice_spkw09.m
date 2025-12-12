function cspice_spkw09( handle, body,  center, frame,  first, ...
                        last,   segid, degree, states, epochs )
   switch nargin
      case 10
         handle = zzmice_int(handle);
         body = zzmice_int(body);
         center = zzmice_int(center);
         frame = zzmice_str(frame);
         first = zzmice_dp(first);
         last = zzmice_dp(last);
         segid = zzmice_str(segid);
         degree = zzmice_int(degree);
         states = zzmice_dp(states);
         epochs = zzmice_dp(epochs);
      otherwise
         error ( ['Usage: cspice_spkw09( handle, body, center, `frame`,' ...
                  ' first, last, `segid`, degree, states(6,n), epochs(n) )'] )
   end
   try
      mice('spkw09_c', handle, body,  center, frame,  first,  ...
                       last,   segid, degree, states, epochs);
   catch spiceerr
      rethrow(spiceerr)
   end