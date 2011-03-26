char const* INPUT_EXAMPLE = "\
\
species : elegans, briggsae, remanei, brenneri\
\
matrix: ase\
A   C   G   T\
-10000 -10000 0 -10000 \
0 -10000 -10000 -10000 \
0 -10000 -10000 -10000 \
-67 -10000 -73 -10000 \
-10000 0 -10000 -10000 \
-10000 -36 -122 -10000 \
-154 -223 -167 -73 \
-167 -201 -182 -67 \
-98 -223 -182 -105 \
-154 -154 -154 -105 \
-154 -223 -154 -78 \
-84 -292 -167 -113\
strand: POS\
total_hits: 500000\
min_score: -1000\
cluster_group: ase\
\
\
consensus_group: urna\
1 GRYNNTCGA\
2 GRYNNTGGA\
cluster_group: ase\
min_score: -1\
\
\
alignment: hlh-2\
ACCTACGCGATC\
CACCGTACGAGT\
ACGACTCGACAG\
AACGATGCATCG\
\
\
\
\
\
position_file: EBRS.gpos \
tree_file: EBRS.itrie \
dna_index_file: EBRS.dnas \
\
cluster: ase\
cluster_type: diversity\
cluster_window: 10000\
cluster_min_requirements: 1\
cluster_requirements_list: [ ase 1 ]\
\
\
hits_file: ase_hits \
cluster_file: ase_clusters \
ones_offset: false \
merge_overlapping_clusters: true\
\
}\
\
\
# make sure you group all things that are used to construct one object\
# in the same rule, so they are parsed together.\
\
";
