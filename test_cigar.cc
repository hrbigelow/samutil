#include "cigar_ops.h"

#include <map>

int main(int argc, char ** argv)
{
    char * genome_to_tx = argv[1];
    char * genome_to_rd = argv[2];

    Cigar::CIGAR_VEC genome_to_tx_cigar = Cigar::FromString(genome_to_tx, 0);
    Cigar::CIGAR_VEC genome_to_rd_cigar = Cigar::FromString(genome_to_rd, 0);

    std::multiset<std::pair<size_t, size_t> > cigar_offsets = 
        Cigar::ComputeOffsets(genome_to_tx_cigar);

    
    Cigar::CIGAR_VEC tx_to_rd_cigar = 
        Cigar::TransitiveMerge(genome_to_tx_cigar, 
                               cigar_offsets,
                               genome_to_rd_cigar, false);

    char tx_to_rd[1000];

    Cigar::ToString(tx_to_rd_cigar.begin(), tx_to_rd_cigar.end(), tx_to_rd);

    fprintf(stdout, "%s\n", tx_to_rd);
    return 0;

}
