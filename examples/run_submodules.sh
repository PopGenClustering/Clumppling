# run clumppling.parseInput to prepare input files
python -m clumppling.parseInput \
-i examples/submodules/input \
-o examples/submodules/output \
-f generalQ 

# run clumppling.alignWithinK on K3 and K5 replicates
python -m clumppling.alignWithinK \
--qfilelist examples/submodules/K3.qfilelist \
-o examples/submodules/K3_aligned.txt

python -m clumppling.alignWithinK \
--qfilelist examples/submodules/K5.qfilelist \
-o examples/submodules/K5_aligned.txt

# run clumppling.detectMode on K3 and K5 modes
python -m clumppling.detectMode \
--align_res examples/submodules/K3_aligned.txt \
--qfilelist examples/submodules/K3.qfilelist \
--qnamelist examples/submodules/K3.qnamelist \
-o examples/submodules/K3_modes

python -m clumppling.detectMode \
--align_res examples/submodules/K5_aligned.txt \
--qfilelist examples/submodules/K5.qfilelist \
--qnamelist examples/submodules/K5.qnamelist \
-o examples/submodules/K5_modes

# prepare mode files
for K in 3 5; do
    SRC=examples/submodules/K${K}_modes
    DST=examples/submodules/K3K5_modes
    mkdir -p $DST
    for f in "$SRC"/*.Q; do
        if [ -f "$f" ]; then
            cp "$f" "$DST/K${K}$(basename "$f")"
        fi
    done
done
# generate list of mode files and names
ls examples/submodules/K3K5_modes/*_rep.Q > examples/submodules/K3K5_modes/K3K5_modes.qfilelist
for f in examples/submodules/K3K5_modes/*_rep.Q; do [ -f "$f" ] && basename "$f" | sed 's/\_rep.Q$//' >> examples/submodules/K3K5_modes/K3K5_modes.qnamelist; done

# run clumppling.alignAcrossK
python -m clumppling.alignAcrossK \
--qfilelist examples/submodules/K3K5_modes/K3K5_modes.qfilelist \
--qnamelist examples/submodules/K3K5_modes/K3K5_modes.qnamelist \
-o examples/submodules/K3K5_acrossK_output