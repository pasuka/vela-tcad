# No-edit Sentaurus Data Explorer round-trip used by the M62 writer gate.
set inp [lindex $cmd_args(rest) 0]
set out [lindex $cmd_args(rest) 1]
if {$inp eq "" || $out eq ""} {
    error "usage: tdx -tcl simplemos_m62_tdx_roundtrip.tcl input.tdr output.tdr"
}
TdrFileOpen $inp
set result [TdrFileSave $inp $out]
TdrFileClose $inp
puts "M62_ROUNDTRIP_OK|$inp|$out|$result"
