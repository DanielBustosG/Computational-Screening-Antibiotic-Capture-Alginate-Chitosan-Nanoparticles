# ==============================================================================
# Clasifica automaticamente cada molecula de antibiotico, en cada frame, como:
#   FREE      = libre en el solvente
#   ADSORBED  = adsorbida/contactando la region superficial
#   ABSORBED  = absorbida en la region interna del polimero
#
# Base del analisis:
#   1) COM dinamico de la nanoparticula (GUL, MAN, CHI).
#   2) Radio dinamico P95 de los atomos pesados del polimero en cada frame.
#   3) Radio interno dinamico = (45/55) * radio P95.
#   4) COM dinamico de cada molecula completa de antibiotico.
#   5) Contacto real por distancia minima entre atomos pesados del antibiotico
#      y atomos pesados del polimero.
#   6) Fraccion de atomos pesados del antibiotico dentro del radio externo/core.
#
# Uso:
#   1. Cargar en VMD la estructura con conectividad y la trayectoria.
#   2. Abrir Tk Console.
#   3. Ejecutar:
#          source classify_antibiotic_location_dynamic.tcl
#
# IMPORTANTE:
#   - La clasificacion usa la molecula completa, no el atomo central auxiliar.
#   - El atomo central se selecciona solo como marcador para plugins posteriores.
#   - El script genera un CSV independiente por antibiotico, para evitar mezclar
#     antibioticos en una misma tabla.
# ============================================================================== 

# ==============================================================================
# 1. PARAMETROS EDITABLES
# ==============================================================================

# Resnames que definen la nanoparticula polimerica.
set polymer_resnames {GUL MAN CHI}

# Los seis antibioticos posibles. El script revisa todos y omite los ausentes.
set antibiotic_resnames {ATC MTC CPZ VGM VMC RFP}

# Numero esperado de copias por tipo. Solo genera advertencia si no coincide.
set expected_molecules_per_antibiotic 8

# Percentil radial que define el radio dinamico externo del polimero.
set polymer_radius_percentile 0.95

# Proporcion geometrica tomada de las regiones iniciales 45 A / 55 A.
set core_fraction [expr {45.0 / 55.0}]

# Distancia maxima para considerar contacto directo entre atomos pesados.
set contact_cutoff_A 3.5

# Cortes adaptativos usados para encontrar la distancia minima real.
# El algoritmo se detiene en el primer corte que contiene al menos un par.
set min_distance_search_cutoffs_A {3.5 5.0 8.0 12.0 20.0 35.0}

# Umbral de control de calidad para marcar alta anisotropia/asfericidad.
# No cambia la categoria; solo genera una advertencia en las salidas.
set shape_anisotropy_warning_threshold 0.30

# Tratamiento de condiciones periodicas.
# 1 = intenta reconstruir y centrar toda la trayectoria con pbctools.
# 0 = usa las coordenadas tal como estan cargadas.
set use_pbc 1

# Tiempo de simulacion.
# Si frame_dt_ns > 0, tiempo = frame * frame_dt_ns.
# Si frame_dt_ns <= 0, se distribuye total_time_ns entre 0 y num_frames-1.
set frame_dt_ns -1.0
set total_time_ns 100.0

# Analizar todos los frames. Cambiar frame_stride para saltar frames.
set first_frame 0
set last_frame -1       ;# -1 = ultimo frame disponible
set frame_stride 1

# Carpeta de resultados.
set output_dir "Antibiotic_location_output"

# ==============================================================================
# 2. FUNCIONES MATEMATICAS Y AUXILIARES
# ==============================================================================

proc distance3 {a b} {
    set dx [expr {[lindex $a 0] - [lindex $b 0]}]
    set dy [expr {[lindex $a 1] - [lindex $b 1]}]
    set dz [expr {[lindex $a 2] - [lindex $b 2]}]
    return [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]
}

proc percentile_nearest_rank {values p} {
    set n [llength $values]
    if {$n == 0} {
        error "No hay valores para calcular el percentil."
    }
    if {$p <= 0.0} { set p 0.0 }
    if {$p >= 1.0} { set p 1.0 }
    set sorted [lsort -real $values]
    set idx [expr {int(ceil($p * $n)) - 1}]
    if {$idx < 0} { set idx 0 }
    if {$idx >= $n} { set idx [expr {$n - 1}] }
    return [lindex $sorted $idx]
}

proc sum_list {values} {
    set total 0.0
    foreach v $values {
        set total [expr {$total + double($v)}]
    }
    return $total
}

proc safe_center_of_mass {sel warning_handle label} {
    set masses [$sel get mass]
    set mass_sum [sum_list $masses]
    if {$mass_sum > 0.0} {
        return [measure center $sel weight mass]
    }
    puts $warning_handle "WARNING: $label no posee masas validas; se usa centro geometrico."
    return [measure center $sel]
}

proc is_hydrogen_atom {element name mass} {
    set e [string toupper [string trim $element]]
    set n [string toupper [string trim $name]]
    if {$e eq "H"} { return 1 }
    if {[string match "H*" $n] && $mass < 2.5} { return 1 }
    return 0
}

proc heavy_indices_from_selection {sel} {
    set indices  [$sel get index]
    set elements [$sel get element]
    set names    [$sel get name]
    set masses   [$sel get mass]
    set heavy {}
    for {set i 0} {$i < [llength $indices]} {incr i} {
        set idx  [lindex $indices $i]
        set elem [lindex $elements $i]
        set name [lindex $names $i]
        set mass [lindex $masses $i]
        if {![is_hydrogen_atom $elem $name $mass]} {
            lappend heavy $idx
        }
    }
    return $heavy
}

proc selection_from_indices {mol indices frame} {
    if {[llength $indices] == 0} {
        error "No se puede crear una seleccion sin indices."
    }

    set sel [atomselect $mol "index [join $indices { }]" frame $frame]

    # La seleccion debe sobrevivir al retorno de esta funcion auxiliar, pero no
    # necesita convertirse en un objeto global de toda la sesion de VMD.
    # "uplevel 1" transfiere la seleccion exactamente al procedimiento que llamo
    # esta funcion:
    #   - si la llama el procedimiento principal, permanece durante todo el loop;
    #   - si la llama adaptive_minimum_distance, permanece solo durante esa
    #     busqueda y se elimina alli mismo.
    # Esto evita tanto "invalid command name atomselectNNN" como el error de
    # bookkeeping interno "can't unset upproc_var_atomselectNNN".
    $sel uplevel 1

    return $sel
}

proc unique_sorted_integers {values} {
    set filtered {}
    foreach v $values {
        if {[string is integer -strict $v] && $v >= 0} {
            lappend filtered $v
        }
    }
    if {[llength $filtered] == 0} { return {} }
    return [lsort -integer -unique $filtered]
}

# Identifica copias moleculares. Primero intenta fragment; si la conectividad no
# separa las copias, usa una llave segname|chain|resid como metodo alternativo.
proc identify_molecule_groups {mol resname expected_count warning_handle} {
    set allsel [atomselect $mol "resname $resname" frame 0]
    if {[$allsel num] == 0} {
        $allsel delete
        return {}
    }

    # Candidatos basados en conectividad (fragment).
    set fragment_groups {}
    set fragments [unique_sorted_integers [$allsel get fragment]]
    if {[llength $fragments] > 1} {
        foreach frag $fragments {
            set s [atomselect $mol "resname $resname and fragment $frag" frame 0]
            set idx [$s get index]
            if {[llength $idx] > 0} {
                lappend fragment_groups [dict create method fragment key $frag indices $idx]
            }
            $s delete
        }
    }

    # Candidatos basados en el indice interno unico de residuo de VMD.
    # Esto es mas robusto que resid cuando varias copias reutilizan el mismo numero.
    set residue_groups {}
    array unset bucket
    set indices [$allsel get index]
    set residues [$allsel get residue]
    for {set i 0} {$i < [llength $indices]} {incr i} {
        set key [lindex $residues $i]
        lappend bucket($key) [lindex $indices $i]
    }
    foreach key [lsort -integer [array names bucket]] {
        lappend residue_groups [dict create method residue_index key $key indices $bucket($key)]
    }

    $allsel delete

    set nf [llength $fragment_groups]
    set nr [llength $residue_groups]

    # Preferir el metodo que coincide con el numero esperado. Si ambos coinciden,
    # fragment conserva la identidad molecular topologica y tiene prioridad.
    if {$nf == $expected_count} {
        return $fragment_groups
    }
    if {$nr == $expected_count} {
        puts $warning_handle "WARNING: $resname se separo mediante residue_index porque fragment detecto $nf grupos."
        return $residue_groups
    }

    # Si ninguno coincide, usar el candidato no trivial mas cercano al esperado.
    if {$nf > 1 && $nr > 1} {
        set diff_f [expr {abs($nf - $expected_count)}]
        set diff_r [expr {abs($nr - $expected_count)}]
        if {$diff_f <= $diff_r} {
            puts $warning_handle "WARNING: $resname usa fragment con $nf grupos; residue_index detecto $nr."
            return $fragment_groups
        } else {
            puts $warning_handle "WARNING: $resname usa residue_index con $nr grupos; fragment detecto $nf."
            return $residue_groups
        }
    }
    if {$nf > 1} { return $fragment_groups }
    if {$nr > 1} { return $residue_groups }

    puts $warning_handle "ERROR: $resname esta presente, pero VMD no pudo separar sus copias moleculares. Revise bonds/fragment y el campo residue. Se omite este antibiotico."
    return {}
}

# Selecciona una sola vez el atomo pesado mas cercano al COM del frame 0.
proc choose_anchor_atom {all_sel heavy_sel warning_handle molecule_id} {
    set com [safe_center_of_mass $all_sel $warning_handle "$molecule_id (anchor)"]
    set indices [$heavy_sel get index]
    set coords  [$heavy_sel get {x y z}]
    set names   [$heavy_sel get name]
    set elements [$heavy_sel get element]

    if {[llength $indices] == 0} {
        error "$molecule_id no tiene atomos pesados para escoger anchor."
    }

    set best_index [lindex $indices 0]
    set best_name [lindex $names 0]
    set best_element [lindex $elements 0]
    set best_distance [distance3 [lindex $coords 0] $com]

    for {set i 1} {$i < [llength $indices]} {incr i} {
        set d [distance3 [lindex $coords $i] $com]
        if {$d < $best_distance} {
            set best_distance $d
            set best_index [lindex $indices $i]
            set best_name [lindex $names $i]
            set best_element [lindex $elements $i]
        }
    }

    return [dict create index $best_index name $best_name element $best_element distance_to_com_A $best_distance]
}

# Calcula radio P95, radio de giro y anisotropia de forma a partir del tensor de
# giro mass-weighted. kappa2 ~ 0: aproximadamente esferica; kappa2 -> 1: alargada.
proc polymer_geometry {heavy_sel np_com radius_percentile warning_handle} {
    set coords [$heavy_sel get {x y z}]
    set masses [$heavy_sel get mass]
    set n [llength $coords]
    if {$n == 0} { error "La seleccion pesada del polimero esta vacia." }

    set mass_sum [sum_list $masses]
    if {$mass_sum <= 0.0} {
        puts $warning_handle "WARNING: masas del polimero no validas para tensor de giro; se usan pesos unitarios."
        set masses [lrepeat $n 1.0]
        set mass_sum [expr {double($n)}]
    }

    set distances {}
    set Sxx 0.0; set Syy 0.0; set Szz 0.0
    set Sxy 0.0; set Sxz 0.0; set Syz 0.0

    for {set i 0} {$i < $n} {incr i} {
        set c [lindex $coords $i]
        set m [lindex $masses $i]
        set dx [expr {[lindex $c 0] - [lindex $np_com 0]}]
        set dy [expr {[lindex $c 1] - [lindex $np_com 1]}]
        set dz [expr {[lindex $c 2] - [lindex $np_com 2]}]
        set r [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]
        lappend distances $r

        set Sxx [expr {$Sxx + $m*$dx*$dx}]
        set Syy [expr {$Syy + $m*$dy*$dy}]
        set Szz [expr {$Szz + $m*$dz*$dz}]
        set Sxy [expr {$Sxy + $m*$dx*$dy}]
        set Sxz [expr {$Sxz + $m*$dx*$dz}]
        set Syz [expr {$Syz + $m*$dy*$dz}]
    }

    foreach var {Sxx Syy Szz Sxy Sxz Syz} {
        set $var [expr {[set $var] / $mass_sum}]
    }

    set trace [expr {$Sxx + $Syy + $Szz}]
    set rg [expr {$trace > 0.0 ? sqrt($trace) : 0.0}]
    set I2 [expr {$Sxx*$Syy + $Syy*$Szz + $Szz*$Sxx - $Sxy*$Sxy - $Sxz*$Sxz - $Syz*$Syz}]

    if {$trace > 0.0} {
        set kappa2 [expr {1.0 - 3.0*$I2/($trace*$trace)}]
    } else {
        set kappa2 0.0
    }
    if {$kappa2 < 0.0 && $kappa2 > -1.0e-10} { set kappa2 0.0 }
    if {$kappa2 > 1.0 && $kappa2 < 1.0 + 1.0e-10} { set kappa2 1.0 }

    set radius [percentile_nearest_rank $distances $radius_percentile]
    return [dict create radius_A $radius rg_A $rg shape_anisotropy $kappa2]
}

proc ligand_radial_metrics {heavy_sel ligand_com np_com np_radius core_radius radius_percentile} {
    set coords [$heavy_sel get {x y z}]
    set n [llength $coords]
    if {$n == 0} { error "Ligando sin atomos pesados." }

    set d_to_lig_com {}
    set inside_np 0
    set inside_core 0

    foreach c $coords {
        lappend d_to_lig_com [distance3 $c $ligand_com]
        set dnp [distance3 $c $np_com]
        if {$dnp <= $np_radius} { incr inside_np }
        if {$dnp <= $core_radius} { incr inside_core }
    }

    set ligand_radius [percentile_nearest_rank $d_to_lig_com $radius_percentile]
    set frac_np [expr {double($inside_np) / double($n)}]
    set frac_core [expr {double($inside_core) / double($n)}]

    return [dict create ligand_radius_A $ligand_radius fraction_inside_np $frac_np fraction_inside_core $frac_core]
}

# Busca la distancia minima real usando measure contacts con cortes crecientes.
# El primer corte que encuentra pares contiene necesariamente el par mas cercano.
proc adaptive_minimum_distance {mol frame lig_heavy_sel polymer_heavy_sel cutoffs} {
    foreach cutoff $cutoffs {
        set contacts [measure contacts $cutoff $lig_heavy_sel $polymer_heavy_sel]
        set list_lig [lindex $contacts 0]
        set list_pol [lindex $contacts 1]
        set npairs [llength $list_lig]

        if {$npairs > 0} {
            # Solo recuperar coordenadas de los atomos que participan en los pares
            # encontrados. Esto evita copiar todos los atomos del polimero para
            # cada antibiotico y mejora mucho el rendimiento.
            set unique_lig [lsort -integer -unique $list_lig]
            set unique_pol [lsort -integer -unique $list_pol]
            set pair_lig_sel [selection_from_indices $mol $unique_lig $frame]
            set pair_pol_sel [selection_from_indices $mol $unique_pol $frame]

            array unset lig_coord
            array unset poly_coord
            set li [$pair_lig_sel get index]
            set lc [$pair_lig_sel get {x y z}]
            set pi [$pair_pol_sel get index]
            set pc [$pair_pol_sel get {x y z}]

            for {set i 0} {$i < [llength $li]} {incr i} {
                set lig_coord([lindex $li $i]) [lindex $lc $i]
            }
            for {set i 0} {$i < [llength $pi]} {incr i} {
                set poly_coord([lindex $pi $i]) [lindex $pc $i]
            }

            set min_d 1.0e30
            for {set i 0} {$i < $npairs} {incr i} {
                set il [lindex $list_lig $i]
                set ip [lindex $list_pol $i]
                set d [distance3 $lig_coord($il) $poly_coord($ip)]
                if {$d < $min_d} { set min_d $d }
            }

            $pair_lig_sel delete
            $pair_pol_sel delete

            return [dict create min_distance_A $min_d search_cutoff_A $cutoff pair_count $npairs found 1]
        }
    }

    # Si no hay pares dentro del maximo buscado, la molecula esta al menos mas
    # lejos que ese corte. La distancia exacta se deja como NA para no crear una
    # busqueda combinatoria innecesaria en moleculas claramente libres.
    return [dict create min_distance_A NA search_cutoff_A [lindex $cutoffs end] pair_count 0 found 0]
}

proc calculate_time_ns {frame first_frame last_frame frame_dt_ns total_time_ns} {
    if {$frame_dt_ns > 0.0} {
        return [expr {($frame - $first_frame) * $frame_dt_ns}]
    }
    if {$last_frame <= $first_frame} { return 0.0 }
    return [expr {double($frame - $first_frame) * $total_time_ns / double($last_frame - $first_frame)}]
}

proc classify_state {dcom core_radius min_distance found_min frac_inside_np contact_cutoff} {
    if {$dcom <= $core_radius} {
        return "ABSORBED"
    }
    set direct_contact 0
    if {$found_min && $min_distance <= $contact_cutoff} {
        set direct_contact 1
    }
    if {$direct_contact || $frac_inside_np > 0.0} {
        return "ADSORBED"
    }
    return "FREE"
}

proc apply_pbc_preprocessing {mol first last polymer_selection warning_handle} {
    if {[catch {package require pbctools} err]} {
        puts $warning_handle "WARNING: pbctools no esta disponible; no se aplico correccion PBC. Detalle: $err"
        return 0
    }

    set cell [molinfo $mol get {a b c}]
    if {[lindex $cell 0] <= 0.0 || [lindex $cell 1] <= 0.0 || [lindex $cell 2] <= 0.0} {
        puts $warning_handle "WARNING: dimensiones periodicas a/b/c no validas; no se aplico correccion PBC."
        return 0
    }

    puts "Applying PBC preprocessing to frames $first:$last ..."
    if {[catch {
        pbc unwrap -molid $mol -first $first -last $last
        pbc wrap -molid $mol -first $first -last $last \
            -compound fragment -center com -centersel $polymer_selection
    } err]} {
        puts $warning_handle "WARNING: fallo la correccion PBC; se continuara con coordenadas cargadas. Detalle: $err"
        return 0
    }
    puts "PBC preprocessing completed."
    return 1
}

# ==============================================================================
# 3. PROGRAMA PRINCIPAL
# ==============================================================================

proc run_antibiotic_location_analysis {} {
    global polymer_resnames antibiotic_resnames expected_molecules_per_antibiotic
    global polymer_radius_percentile core_fraction contact_cutoff_A
    global min_distance_search_cutoffs_A shape_anisotropy_warning_threshold
    global use_pbc frame_dt_ns total_time_ns first_frame last_frame frame_stride
    global output_dir

    if {[molinfo top] == -1} {
        puts "ERROR: no hay una molecula/trajectory cargada en VMD."
        return
    }

    set mol [molinfo top]
    set num_frames_total [molinfo $mol get numframes]
    if {$num_frames_total < 1} {
        puts "ERROR: la molecula no contiene frames."
        return
    }

    set local_last $last_frame
    if {$local_last < 0 || $local_last >= $num_frames_total} {
        set local_last [expr {$num_frames_total - 1}]
    }
    if {$first_frame < 0 || $first_frame > $local_last} {
        puts "ERROR: rango de frames invalido: first=$first_frame last=$local_last"
        return
    }
    if {$frame_stride < 1} {
        puts "ERROR: frame_stride debe ser >= 1."
        return
    }

    if {![file exists $output_dir]} {
        file mkdir $output_dir
    }

    set warning_path [file join $output_dir "analysis_warnings.log"]
    set warning_handle [open $warning_path w]
    puts $warning_handle "Antibiotic location analysis warnings"
    puts $warning_handle "Molecule id: $mol"

    set polymer_selection "resname [join $polymer_resnames { }]"
    set polymer_all [atomselect $mol $polymer_selection frame $first_frame]
    if {[$polymer_all num] == 0} {
        puts "ERROR: no se encontraron atomos para: $polymer_selection"
        puts $warning_handle "ERROR: seleccion polimerica vacia: $polymer_selection"
        close $warning_handle
        return
    }

    set polymer_heavy_indices [heavy_indices_from_selection $polymer_all]
    if {[llength $polymer_heavy_indices] == 0} {
        puts "ERROR: no se encontraron atomos pesados en la nanoparticula."
        puts $warning_handle "ERROR: sin atomos pesados en nanoparticula."
        $polymer_all delete
        close $warning_handle
        return
    }
    set polymer_heavy [selection_from_indices $mol $polymer_heavy_indices $first_frame]

    puts "============================================================"
    puts "Dynamic antibiotic location analysis"
    puts "Molecule id: $mol"
    puts "Frames loaded: $num_frames_total"
    puts "Frames analyzed: $first_frame to $local_last, stride $frame_stride"
    puts "Polymer selection: $polymer_selection"
    puts "Polymer atoms: [$polymer_all num]"
    puts "Polymer heavy atoms: [$polymer_heavy num]"
    puts "============================================================"

    if {$use_pbc} {
        apply_pbc_preprocessing $mol $first_frame $local_last $polymer_selection $warning_handle
    }

    # Detectar y preparar antibioticos.
    set ligand_data [dict create]
    set detected_handle [open [file join $output_dir "antibiotics_detected.txt"] w]
    set anchor_handle [open [file join $output_dir "central_atom_indices.txt"] w]

    puts $detected_handle "Molecule_id: $mol"
    puts $detected_handle "Frames_analyzed: $first_frame-$local_last stride=$frame_stride"
    puts $detected_handle "Polymer_selection: $polymer_selection"
    puts $detected_handle "Radius_definition: heavy-atom P[expr {int($polymer_radius_percentile*100)}]"
    puts $detected_handle "Core_fraction: $core_fraction"
    puts $detected_handle "Contact_cutoff_A: $contact_cutoff_A"
    puts $detected_handle ""

    puts $anchor_handle "One fixed representative heavy atom per antibiotic molecule"
    puts $anchor_handle "These atoms are auxiliary markers only; classification uses the complete molecule."
    puts $anchor_handle ""

    set all_anchor_indices {}
    set detected_antibiotics {}

    foreach rn $antibiotic_resnames {
        set testsel [atomselect $mol "resname $rn" frame $first_frame]
        set natoms [$testsel num]
        $testsel delete

        if {$natoms == 0} {
            puts "$rn: absent; skipped."
            puts $detected_handle "$rn: absent"
            continue
        }

        set raw_groups [identify_molecule_groups $mol $rn $expected_molecules_per_antibiotic $warning_handle]
        set ngroups [llength $raw_groups]
        if {$ngroups == 0} {
            puts "$rn: present, but molecule copies could not be separated; skipped."
            puts $detected_handle "$rn: present, separation_failed"
            continue
        }

        if {$ngroups != $expected_molecules_per_antibiotic} {
            puts $warning_handle "WARNING: $rn expected=$expected_molecules_per_antibiotic detected=$ngroups. Analysis continues."
        }

        set prepared_groups {}
        set group_number 0
        foreach raw $raw_groups {
            incr group_number
            set indices [dict get $raw indices]
            set all_sel [selection_from_indices $mol $indices $first_frame]
            set heavy_indices [heavy_indices_from_selection $all_sel]
            if {[llength $heavy_indices] == 0} {
                puts $warning_handle "WARNING: $rn group $group_number has no heavy atoms; omitted."
                $all_sel delete
                continue
            }
            set heavy_sel [selection_from_indices $mol $heavy_indices $first_frame]
            set molecule_id [format "%s_%02d" $rn $group_number]
            set anchor [choose_anchor_atom $all_sel $heavy_sel $warning_handle $molecule_id]
            set anchor_index [dict get $anchor index]
            lappend all_anchor_indices $anchor_index

            puts $anchor_handle [format "%s,resname=%s,method=%s,key=%s,anchor_index=%s,anchor_name=%s,anchor_element=%s,distance_to_COM_A=%.6f" \
                $molecule_id $rn [dict get $raw method] [dict get $raw key] \
                $anchor_index [dict get $anchor name] [dict get $anchor element] \
                [dict get $anchor distance_to_com_A]]

            lappend prepared_groups [dict create \
                molecule_id $molecule_id \
                method [dict get $raw method] \
                key [dict get $raw key] \
                all_sel $all_sel \
                heavy_sel $heavy_sel \
                anchor_index $anchor_index \
                anchor_name [dict get $anchor name]]
        }

        if {[llength $prepared_groups] == 0} {
            puts "$rn: no valid molecular groups after validation; skipped."
            puts $detected_handle "$rn: present, no_valid_groups"
            continue
        }

        dict set ligand_data $rn $prepared_groups
        lappend detected_antibiotics $rn
        puts "$rn: present, [llength $prepared_groups] molecules detected."
        puts $detected_handle "$rn: present, molecules=[llength $prepared_groups]"
    }

    puts $anchor_handle ""
    puts $anchor_handle "VMD_selection: index [join $all_anchor_indices { }]"
    close $anchor_handle
    close $detected_handle

    if {[llength $detected_antibiotics] == 0} {
        puts "ERROR: ninguno de los antibioticos esperados pudo analizarse."
        puts $warning_handle "ERROR: no analyzable antibiotics detected."
        $polymer_all delete
        $polymer_heavy delete
        close $warning_handle
        return
    }

    # Abrir archivos independientes por antibiotico.
    array unset positions_handle
    array unset counts_handle
    foreach rn $detected_antibiotics {
        set positions_handle($rn) [open [file join $output_dir "${rn}_positions.csv"] w]
        puts $positions_handle($rn) "Frame,Time_ns,Molecule_ID,Grouping_method,Grouping_key,Anchor_index,Anchor_name,Ligand_COM_x_A,Ligand_COM_y_A,Ligand_COM_z_A,NP_COM_x_A,NP_COM_y_A,NP_COM_z_A,COM_distance_A,NP_radius_P95_A,NP_diameter_A,Core_radius_A,Normalized_distance,Penetration_A,Ligand_radius_P95_A,Radial_surface_gap_A,Min_polymer_distance_A,Min_distance_search_cutoff_A,Contact_pair_count,Direct_contact,Fraction_inside_NP,Fraction_inside_core,State,NP_Rg_A,Shape_anisotropy,Shape_warning,Radial_compaction,Approx_volume_ratio"

        set counts_handle($rn) [open [file join $output_dir "${rn}_state_counts.csv"] w]
        puts $counts_handle($rn) "Frame,Time_ns,Free,Adsorbed,Absorbed,Total"
    }

    set geometry_handle [open [file join $output_dir "nanoparticle_geometry.csv"] w]
    puts $geometry_handle "Frame,Time_ns,NP_COM_x_A,NP_COM_y_A,NP_COM_z_A,NP_radius_P95_A,NP_diameter_A,Core_radius_A,NP_Rg_A,Shape_anisotropy,Shape_warning,Radial_compaction,Approx_volume_ratio"

    # Acumuladores para resumen por molecula.
    array unset state_free
    array unset state_adsorbed
    array unset state_absorbed
    array unset first_adsorption
    array unset first_absorption
    array unset total_observations

    set analyzed_frames {}
    for {set f $first_frame} {$f <= $local_last} {incr f $frame_stride} {
        lappend analyzed_frames $f
    }
    set n_analyzed_frames [llength $analyzed_frames]

    set initial_radius -1.0
    set processed 0

    foreach f $analyzed_frames {
        incr processed
        set time_ns [calculate_time_ns $f $first_frame $local_last $frame_dt_ns $total_time_ns]

        $polymer_all frame $f
        $polymer_all update
        $polymer_heavy frame $f
        $polymer_heavy update

        set np_com [safe_center_of_mass $polymer_all $warning_handle "polymer frame $f"]
        set geom [polymer_geometry $polymer_heavy $np_com $polymer_radius_percentile $warning_handle]
        set np_radius [dict get $geom radius_A]
        set np_rg [dict get $geom rg_A]
        set shape [dict get $geom shape_anisotropy]
        set np_diameter [expr {2.0 * $np_radius}]
        set core_radius [expr {$core_fraction * $np_radius}]

        if {$initial_radius < 0.0} { set initial_radius $np_radius }
        if {$initial_radius > 0.0} {
            set radial_compaction [expr {$np_radius / $initial_radius}]
            set volume_ratio [expr {pow($radial_compaction, 3.0)}]
        } else {
            set radial_compaction 1.0
            set volume_ratio 1.0
        }

        if {$shape >= $shape_anisotropy_warning_threshold} {
            set shape_warning "HIGH_ASPHERICITY"
            puts $warning_handle [format "WARNING: frame=%d shape_anisotropy=%.6f exceeds threshold=%.3f" $f $shape $shape_anisotropy_warning_threshold]
        } else {
            set shape_warning "OK"
        }

        puts $geometry_handle [format "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%.6f,%.6f" \
            $f $time_ns [lindex $np_com 0] [lindex $np_com 1] [lindex $np_com 2] \
            $np_radius $np_diameter $core_radius $np_rg $shape $shape_warning \
            $radial_compaction $volume_ratio]

        foreach rn $detected_antibiotics {
            set count_free 0
            set count_adsorbed 0
            set count_absorbed 0

            foreach group [dict get $ligand_data $rn] {
                set all_sel [dict get $group all_sel]
                set heavy_sel [dict get $group heavy_sel]
                $all_sel frame $f
                $all_sel update
                $heavy_sel frame $f
                $heavy_sel update

                set molecule_id [dict get $group molecule_id]
                set ligand_com [safe_center_of_mass $all_sel $warning_handle "$molecule_id frame $f"]
                set dcom [distance3 $ligand_com $np_com]
                set normalized_distance [expr {$np_radius > 0.0 ? $dcom / $np_radius : 0.0}]
                set penetration [expr {$np_radius - $dcom}]

                set radial [ligand_radial_metrics $heavy_sel $ligand_com $np_com $np_radius $core_radius $polymer_radius_percentile]
                set ligand_radius [dict get $radial ligand_radius_A]
                set frac_np [dict get $radial fraction_inside_np]
                set frac_core [dict get $radial fraction_inside_core]

                set mindata [adaptive_minimum_distance $mol $f $heavy_sel $polymer_heavy $min_distance_search_cutoffs_A]
                set found_min [dict get $mindata found]
                set min_distance [dict get $mindata min_distance_A]
                set search_cutoff [dict get $mindata search_cutoff_A]
                set pair_count [dict get $mindata pair_count]

                set radial_surface_gap [expr {$dcom - $np_radius - $ligand_radius}]
                set direct_contact [expr {$found_min && $min_distance <= $contact_cutoff_A ? 1 : 0}]
                set state [classify_state $dcom $core_radius $min_distance $found_min $frac_np $contact_cutoff_A]

                switch -- $state {
                    FREE { incr count_free }
                    ADSORBED { incr count_adsorbed }
                    ABSORBED { incr count_absorbed }
                }

                set key "$rn|$molecule_id"
                if {![info exists state_free($key)]} {
                    set state_free($key) 0
                    set state_adsorbed($key) 0
                    set state_absorbed($key) 0
                    set first_adsorption($key) "NA"
                    set first_absorption($key) "NA"
                    set total_observations($key) 0
                }
                incr total_observations($key)
                switch -- $state {
                    FREE { incr state_free($key) }
                    ADSORBED {
                        incr state_adsorbed($key)
                        if {$first_adsorption($key) eq "NA"} { set first_adsorption($key) $time_ns }
                    }
                    ABSORBED {
                        incr state_absorbed($key)
                        if {$first_absorption($key) eq "NA"} { set first_absorption($key) $time_ns }
                    }
                }

                if {$found_min} {
                    set min_distance_text [format "%.6f" $min_distance]
                } else {
                    set min_distance_text "NA"
                }

                puts $positions_handle($rn) [format "%d,%.6f,%s,%s,%s,%s,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%.6f,%d,%d,%.6f,%.6f,%s,%.6f,%.6f,%s,%.6f,%.6f" \
                    $f $time_ns $molecule_id [dict get $group method] [dict get $group key] \
                    [dict get $group anchor_index] [dict get $group anchor_name] \
                    [lindex $ligand_com 0] [lindex $ligand_com 1] [lindex $ligand_com 2] \
                    [lindex $np_com 0] [lindex $np_com 1] [lindex $np_com 2] \
                    $dcom $np_radius $np_diameter $core_radius $normalized_distance \
                    $penetration $ligand_radius $radial_surface_gap $min_distance_text \
                    $search_cutoff $pair_count $direct_contact $frac_np $frac_core $state \
                    $np_rg $shape $shape_warning $radial_compaction $volume_ratio]
            }

            set total_count [expr {$count_free + $count_adsorbed + $count_absorbed}]
            puts $counts_handle($rn) [format "%d,%.6f,%d,%d,%d,%d" \
                $f $time_ns $count_free $count_adsorbed $count_absorbed $total_count]
        }

        if {$processed == 1 || $processed % 10 == 0 || $processed == $n_analyzed_frames} {
            puts [format "Processed frame %d (%d/%d), time=%.3f ns, R95=%.3f A, kappa2=%.4f" \
                $f $processed $n_analyzed_frames $time_ns $np_radius $shape]
        }
    }

    # Cerrar archivos de trayectorias y conteos.
    foreach rn $detected_antibiotics {
        close $positions_handle($rn)
        close $counts_handle($rn)
    }
    close $geometry_handle

    # Resumen independiente por antibiotico.
    foreach rn $detected_antibiotics {
        set summary_handle [open [file join $output_dir "${rn}_summary.csv"] w]
        puts $summary_handle "Molecule_ID,Free_frames,Adsorbed_frames,Absorbed_frames,Total_frames,Free_percent,Adsorbed_percent,Absorbed_percent,First_adsorption_ns,First_absorption_ns"
        foreach group [dict get $ligand_data $rn] {
            set molecule_id [dict get $group molecule_id]
            set key "$rn|$molecule_id"
            set total $total_observations($key)
            set free_n $state_free($key)
            set ads_n $state_adsorbed($key)
            set abs_n $state_absorbed($key)
            set free_pct [expr {100.0 * $free_n / double($total)}]
            set ads_pct [expr {100.0 * $ads_n / double($total)}]
            set abs_pct [expr {100.0 * $abs_n / double($total)}]
            puts $summary_handle [format "%s,%d,%d,%d,%d,%.6f,%.6f,%.6f,%s,%s" \
                $molecule_id $free_n $ads_n $abs_n $total $free_pct $ads_pct $abs_pct \
                $first_adsorption($key) $first_absorption($key)]
        }
        close $summary_handle
    }

    # Liberar selecciones.
    foreach rn $detected_antibiotics {
        foreach group [dict get $ligand_data $rn] {
            set s1 [dict get $group all_sel]
            set s2 [dict get $group heavy_sel]
            $s1 delete
            $s2 delete
        }
    }
    $polymer_all delete
    $polymer_heavy delete

    puts $warning_handle "Analysis completed."
    close $warning_handle

    puts "============================================================"
    puts "Analysis completed successfully."
    puts "Detected antibiotics: [join $detected_antibiotics { }]"
    puts "Results directory: $output_dir"
    puts "One positions CSV, one state-count CSV, and one summary CSV were"
    puts "generated for each detected antibiotic."
    puts "Use central_atom_indices.txt for the fixed anchor selection."
    puts "============================================================"
}

# Ejecutar automaticamente al hacer source y mostrar una traza util si VMD
# encuentra un error. Esto facilita localizar la linea exacta sin dejar el
# diagnostico reducido al ultimo mensaje de Tcl.
if {[catch {run_antibiotic_location_analysis} analysis_error analysis_options]} {
    puts stderr ""
    puts stderr "============================================================"
    puts stderr "ERROR: antibiotic location analysis stopped."
    puts stderr "Message: $analysis_error"
    if {[dict exists $analysis_options -errorinfo]} {
        puts stderr "Detailed Tcl trace:"
        puts stderr [dict get $analysis_options -errorinfo]
    }
    puts stderr "============================================================"
}
