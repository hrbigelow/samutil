#include <cstdio>
#include <cstdlib>

#include "pileup_tools.h"
#include "base_qual_strand_reader.h"
#include "nucleotide_stats.h"

int bqs_usage()
{
    fprintf(stderr,
            "Usage: dep bqs input.pileup > counts.bqs\n");
    return 1;
}

int main_bqs(int argc, char ** argv)
{

    if (argc != 2)
    {
        return bqs_usage();
    }

    char * pileup_input_file = argv[1];
    size_t const MAX_QUALITY = 255;

    char const* strands = "+-";

    BaseQualStrandReader reader;
    reader.initialize(pileup_input_file);

    double fake_nuc_frequency[] = { 1, 1, 1, 1 };

    //all we need is a name mapping.
    size_t datum_index = 0;

    JPD_DATA fake_jpd;
    for (size_t b = 0; b != 4; ++b)
    {
        for (size_t q = 0; q <= MAX_QUALITY; ++q)
        {
            for (size_t s = 0; s != 2; ++s)
            {
                BaseQualStrandReader::Datum datum = 
                    { Nucleotide::bases_upper[b], q, strands[s], datum_index, 0 };

                fake_jpd.insert(std::make_pair(datum.name(),nuc_frequency(fake_nuc_frequency)));
                ++datum_index;
            }
        }
    }

    NucleotideStats fake_stats(fake_jpd.size());
    fake_stats.initialize(fake_jpd);

    double * all_counts = new double[fake_jpd.size()];
    std::fill(all_counts, all_counts + fake_jpd.size(), 0.0);

    size_t min_quality_score = 0;

    while (reader.more_loci())
    {
        LocusSummary locus = 
            reader.get_next_locus(fake_stats, static_cast<void const*>(& min_quality_score));

        for (size_t raw_index = 0; raw_index != locus.num_distinct_data; ++raw_index)
        {
            size_t data_index = locus.stats_index[raw_index];
            all_counts[data_index] += locus.raw_counts[raw_index];
        }
    }


    size_t max_quality_present = 0;

    std::map<std::string, size_t>::const_iterator datum_iter;

    //find the maximum nonzero quality present
    for (datum_iter = fake_stats.name_mapping.begin();
         datum_iter != fake_stats.name_mapping.end();
         ++datum_iter)
    {
        BaseQualStrandReader::Datum datum = 
            reader.get_datum_from_name((*datum_iter).first);
        size_t datum_index = (*datum_iter).second;

        if (all_counts[datum_index] > 0.0)
        {
            max_quality_present = std::max(max_quality_present, datum.quality);
        }
        
    }

    for (datum_iter = fake_stats.name_mapping.begin();
         datum_iter != fake_stats.name_mapping.end();
         ++datum_iter)
    {
        BaseQualStrandReader::Datum datum = 
            reader.get_datum_from_name((*datum_iter).first);
        size_t datum_index = (*datum_iter).second;

        if (datum.quality <= max_quality_present)
        {
            printf("%c\t%Zu\t%c\t%Zu\n",
                   datum.called_base,
                   datum.quality,
                   datum.strand,
                   static_cast<size_t>(all_counts[datum_index]));
        }
    }
    
    delete all_counts;
    
    return 0;
}

