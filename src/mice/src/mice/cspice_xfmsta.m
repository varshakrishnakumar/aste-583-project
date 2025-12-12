function [ostate] = cspice_xfmsta( istate, icosys, ocosys, body )
   switch nargin
      case 4
         istate = zzmice_dp (istate);
         icosys = zzmice_str(icosys);
         ocosys = zzmice_str(ocosys);
         body   = zzmice_str(body);
      otherwise
         error ( ['Usage: [_ostate(6)_] = cspice_xfmsta( `_istate(6)_`, ' ...
                  '`icosys`, `ocosys`, `body`)'] )
   end
   try
      [ostate] = mice('xfmsta_c', istate, icosys, ocosys, body);
   catch spiceerr
      rethrow(spiceerr)
   end