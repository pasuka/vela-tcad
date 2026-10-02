# Read-only Sentaurus Data Explorer inventory used by SimpleMOS M62.
set inp [lindex $cmd_args(rest) 0]
if {$inp eq ""} {
    error "usage: tdx -tcl simplemos_m62_tdx_inventory.tcl input.tdr"
}
TdrFileOpen $inp
set ng [TdrFileGetNumGeometry $inp]
puts "record|geometry|region|state|dataset|region_name|state_name|dataset_name|quantity|location|structure|value_count"
for {set ig 0} {$ig < $ng} {incr ig} {
    set ns [TdrGeometryGetNumState $inp $ig]
    set nr [TdrGeometryGetNumRegion $inp $ig]
    for {set is 0} {$is < $ns} {incr is} {
        set sname [TdrStateGetName $inp $ig $is]
        for {set ir 0} {$ir < $nr} {incr ir} {
            set rname [TdrRegionGetName $inp $ig $ir]
            set nd [TdrRegionGetNumDataset $inp $ig $ir $is]
            for {set id 0} {$id < $nd} {incr id} {
                set dname [TdrDatasetGetName $inp $ig $ir $is $id]
                set quantity [TdrDatasetGetQuantity $inp $ig $ir $is $id]
                set location [TdrDatasetGetLocation $inp $ig $ir $is $id]
                set structure [TdrDatasetGetStructure $inp $ig $ir $is $id]
                set nv [TdrDatasetGetNumValue $inp $ig $ir $is $id]
                puts "dataset|$ig|$ir|$is|$id|$rname|$sname|$dname|$quantity|$location|$structure|$nv"
            }
        }
    }
}
TdrFileClose $inp
