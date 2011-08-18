#define __STDC_LIMIT_MACROS
#include <pthread.h>
#include <stdint.h>
#include <algorithm>
#include <map>
#include <set>

#include "readsim_aux.h"
#include "cisortho/nested.h"
#include "dep/tools.h"
#include "cisortho/region.h"
#include "sam_buffer.h"
#include "sam_helper.h"
#include "sam_score_aux.h"


/*

projects transcriptome-aligned reads to the genome.

transcript -> genome, transcript -> read :    genome -> read

NEEDED: 

1. tx_to_genome  (projections of transcripts to the genome)
2. tx_name_to_projection.  (tx names pointing to tx_to_genome iterators)

0.  source = transcript, target = genome
1.  parse a samline
2.  retrieve the source_to_target cigar from the tx_to_genome 
3.  merge the source_to_target and source_to_read CIGARs to create a target_to_read CIGAR

If our input is samlines aligned to transcripts, and the transcripts are ordered
by start on the genome, then we just keep track of new transcripts in the SamLine input.

*/
int transcript_to_genome_usage()
{
    fprintf(stderr,
            "\nUsage:\n\n"
            "samutil tx2genome [OPTIONS] transcripts.gtf reads.vs.tx.asort.sam reads.vs.genome.sam\n\n"
            "Options:\n\n"
            "-h     STRING   input genome SAM header file\n"
            "-n     FLAG     If present, use CIGAR 'N' to represent introns.  Otherwise use 'D'\n"
            "-x     FLAG     If present, add Cufflinks XS:A: tag (+/-) denotes source RNA strand.\n"
            "-t     INT      Number of threads to use\n"
            "\n"
            "Only one of each group of identical fragment alignments is output.\n"
            "These arise from congruent subsets of isoforms in transcripts.gtf.\n"
            "SAM Records successfully projected will have the 'XP:A:T' tag added.\n"
            "reads.vs.tx.sam must be sorted by [rname, read_pair_flag, pos].\n"
            "\n"
            "reads.vs.tx.asort.sam must be sorted by alignment position, with contigs (transcript)\n"
            "ordered by projection order according to genomic contig order as specified in\n\n"
            );
    return 1;
}

typedef std::set<SequenceProjection>::const_iterator SP_ITER;
typedef std::map<char const*, SP_ITER, less_char_ptr> PROJ_MAP;
typedef PROJ_MAP::const_iterator NP_ITER;



struct project_aux_data
{
    project_aux_data(char const* _p,     
                     PAIRED_READ_SET::const_iterator _start,
                     PAIRED_READ_SET::const_iterator _stop,
                     SP_ITER _f,
                     SP_ITER _s,
                     PROJ_MAP const* _t,
                     NP_ITER const* _n,
                     bool _i,
                     bool _ac,
                     SamOrder const* _rec,
                     CONTIG_OFFSETS::const_iterator _ctg) :
        prev_rname(_p),
        start(_start),
        stop(_stop),
        first_seq_proj(_f),
        second_seq_proj(_s),
        tx_name_to_projection(_t),
        name_to_proj(_n),
        inserts_are_introns(_i),
        add_cufflinks_xs_tag(_ac),
        record_ordering(_rec),
        contig_iter(_ctg) { }

    project_aux_data() { }

    char const* prev_rname;
    PAIRED_READ_SET::const_iterator start;
    PAIRED_READ_SET::const_iterator stop;
    SP_ITER first_seq_proj; //defines the transcript-to-genome sequence projection
    SP_ITER second_seq_proj;
    PROJ_MAP const* tx_name_to_projection;
    NP_ITER const* name_to_proj; //shortcut for lookup of second projection in PROJ_MAP
    bool inserts_are_introns;
    bool add_cufflinks_xs_tag;
    SamOrder const* record_ordering;
    CONTIG_OFFSETS::const_iterator contig_iter;
};



void * project_transcript_entries_aux(void * data)
{
    project_aux_data & p = *static_cast<project_aux_data *>(data);

    for (PAIRED_READ_SET::const_iterator pit = p.start; pit != p.stop; ++pit)
    {
        SamLine * first = const_cast<SamLine *>((*pit).first);
        SamLine * second = const_cast<SamLine *>((*pit).second);

        assert(strcmp(p.prev_rname, first->rname) == 0);

        //but, their mates may not be
        NP_ITER nit_mate = strcmp(first->rname, first->mate_ref_name()) == 0 ? *p.name_to_proj
            : p.tx_name_to_projection->find(first->mate_ref_name());
                    
        if (nit_mate == p.tx_name_to_projection->end())
        {
            fprintf(stderr, "Error: couldn't find transcript %s in projections "
                    "derived from GTF file\n", first->mate_ref_name());
            exit(1);
        }
        else
        {
            p.second_seq_proj = (*nit_mate).second;
        }
                    
        bool projection_applied = 
            ApplyProjectionToSAM(*p.first_seq_proj, *p.second_seq_proj, "T", 
                                 first, second, 
                                 p.inserts_are_introns,
                                 p.add_cufflinks_xs_tag);
                        
        assert(projection_applied);

        assert(AreMappedMatePairs(*first, *second));

        first->SetFlattenedPosition(p.record_ordering->contig_offsets, &p.contig_iter);
        second->SetFlattenedPosition(p.record_ordering->contig_offsets, &p.contig_iter);
    }
    pthread_exit(0);
}


//call at the end of loading the input_buffer with all alignments to a
//given transcript.

//find projection by name.  assume all unique entry pairs in input
//buffer are on the same transcript, namely 'prev_rname'. Project them
//to the output buffer, and update the lowbound for the output buffer
//to the minimum position on the genome of the transcript being processed.
void ProjectTranscriptEntries(PROJ_MAP const& tx_name_to_projection, 
                              SamBuffer & input_buffer, 
                              SamBuffer & output_buffer,
                              char const* prev_rname,
                              SamLine * low_bound,
                              bool inserts_are_introns,
                              bool add_cufflinks_xs_tag,
                              size_t num_threads,
                              project_aux_data * thread_data_array)
{
    if (input_buffer.unique_entry_pairs.empty())
    {
        return;
    }

    std::pair<SamLine const*, bool> insert_result;

    SP_ITER first_proj;
    SP_ITER second_proj;
    NP_ITER nit = tx_name_to_projection.find(prev_rname);
    if (nit == tx_name_to_projection.end())
    {
        fprintf(stderr, "Error: couldn't find transcript %s in projections "
                "derived from GTF file\n", prev_rname);
        exit(1);
    }
    else
    {
        first_proj = (*nit).second;
    }

    CONTIG_OFFSETS::const_iterator contig_iter = 
        output_buffer.sam_order->contig_offsets.begin();

    //now, 

    size_t small_chunk_size = input_buffer.unique_entry_pairs.size() / num_threads;
    size_t num_plus_one_threads = input_buffer.unique_entry_pairs.size() % num_threads;
    size_t num_nonzero_threads = small_chunk_size == 0 ? num_plus_one_threads : num_threads;

    PAIRED_READ_SET::iterator start = input_buffer.unique_entry_pairs.begin();
    PAIRED_READ_SET::iterator stop = start;

    int * thread_retval = new int[num_threads];
    pthread_t * thread_ids = new pthread_t[num_threads];

    for (size_t c = 0; c != num_nonzero_threads; ++c)
    {
        size_t this_chunk_size = 
            (c < num_plus_one_threads ? small_chunk_size + 1 : small_chunk_size);

        std::advance(stop, this_chunk_size);

        project_aux_data * p = & thread_data_array[c];

        p->prev_rname = prev_rname;
        p->start = start;
        p->stop = stop;
        p->first_seq_proj = first_proj;
        p->second_seq_proj = second_proj;
        p->name_to_proj = &nit;
        p->contig_iter = contig_iter;

        thread_retval[c] = 
            pthread_create(&thread_ids[c], NULL, project_transcript_entries_aux, 
                           static_cast<void *>(p));

        start = stop;
    }

    for (size_t c = 0; c != num_nonzero_threads; ++c)
    {
        pthread_join(thread_ids[c], NULL);
    }

    delete thread_retval;
    delete thread_ids;

    //parallelize this with a fork.  How?
    for (PAIRED_READ_SET::iterator pit = input_buffer.unique_entry_pairs.begin();
         pit != input_buffer.unique_entry_pairs.end(); ++pit)
    {

        insert_result = output_buffer.insert((*pit).first);
        insert_result = output_buffer.insert((*pit).second);

    }
    input_buffer.unique_entry_pairs.clear();
                
    //at this point, the low bound
    char const* target_contig = (*first_proj).target_dna.c_str();
    size_t target_min_pos = (*first_proj).cigar[0].length;

    
    strcpy(low_bound->rname, target_contig);
    low_bound->pos = target_min_pos;

    CONTIG_OFFSETS::const_iterator iter_tmp(thread_data_array[0].contig_iter);
    low_bound->SetFlattenedPosition(output_buffer.sam_order->contig_offsets, &iter_tmp);

}


//check that the genome contig order implied by the combination of
//transcript projections and transcript contig order is the same as
//that given by genome_contig_order
bool CheckProjectionOrder(CONTIG_OFFSETS const& genome_contig_order,
                          CONTIG_OFFSETS const& tx_contig_order,
                          std::set<SequenceProjection> const& tx_to_genome)
{
    //1. Compute a flattened genome coordinate for each transcript
    std::multimap<size_t, char const*> tx_implied_genomic_offsets;
    std::multimap<size_t, char const*>::const_iterator tit;

    std::set<SequenceProjection>::const_iterator pit;
    OFFSETS_ITER oit;
    for (pit = tx_to_genome.begin(); pit != tx_to_genome.end(); ++pit)
    {
        char const* tx_contig = (*pit).source_dna.c_str();
        char const* tx_proj_contig = (*pit).target_dna.c_str();
        oit = genome_contig_order.find(tx_proj_contig);
        if (oit == genome_contig_order.end())
        {
            fprintf(stderr, "Error: CheckProjectionOrder: didn't find genome contig %s "
                    "that was found in transcript-to-genome projection from gtf file\n",
                    tx_proj_contig);
            exit(1);
        }
        size_t genome_contig_offset = (*oit).second;
        size_t tx_local_offset = Cigar::LeftOffset((*pit).cigar, false);

        //assert(tx_implied_genomic_offsets.find() == tx_implied_genomic_offsets.end());
        
        tx_implied_genomic_offsets.insert(std::make_pair(genome_contig_offset + tx_local_offset, tx_contig));
    }

    size_t last_implied_offset = SIZE_MAX;
    size_t last_given_offset = 0;
    char const* last_tx_contig = "";

    //tx implied offsets are partially ordered.  given offsets are
    //fully ordered.  here, we merely check that if a consecutive pair
    //of contigs are implied strictly increasing, they are increasing
    //by given offset as well. 
    for (tit = tx_implied_genomic_offsets.begin(); tit != tx_implied_genomic_offsets.end(); ++tit)
    {
        size_t implied_offset = (*tit).first;
        char const* tx_contig = (*tit).second;
        oit = tx_contig_order.find(tx_contig);
        if (oit == tx_contig_order.end())
        {
            fprintf(stderr, "Error: CheckProjectionOrder: didn't find transcript contig %s "
                    "in given transcript ordering that was found in transcript-to-genome "
                    "projection from gtf file\n",
                    tx_contig);
            exit(1);
        }
        size_t given_offset = (*oit).second;
        if (last_implied_offset < implied_offset
            && ! (last_given_offset < given_offset))
        {
            fprintf(stderr, "Error: CheckProjectionOrder: transcript contig ordering implied from "
                    "gtf file and given genome contig ordering differs from "
                    "given transcript ordering. Please re-order transcript contigs in SAM header\n"
                    "Last contig: %s (%Zu implied offset) (%Zu given offset)\n"
                    "Current ctg: %s (%Zu implied offset) (%Zu given offset)\n",
                    last_tx_contig, last_implied_offset, last_given_offset, 
                    tx_contig, implied_offset, given_offset);

            return false;
        }
        else
        {
            last_implied_offset = implied_offset;
            last_given_offset = given_offset;
            last_tx_contig = tx_contig;
        }
    }
    return true;
}




int main_transcript_to_genome(int argc, char ** argv)
{
    char c;

    char const* input_genome_sam_header_file = "";
    bool inserts_are_introns = false;
    bool add_cufflinks_xs_tag = false;
    size_t num_threads = 1;

    while ((c = getopt(argc, argv, "h:nxt:")) >= 0)
    {
        switch(c)
        {
        case 'h': input_genome_sam_header_file = optarg; break;
        case 'n': inserts_are_introns = true; break;
        case 'x': add_cufflinks_xs_tag = true; break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        default: return transcript_to_genome_usage(); break;
        }
    }

    //printf("optind: %i\n", optind);
    if (argc != 3 + optind)
    {
        return transcript_to_genome_usage();
    }

    char * gtf_file = argv[optind];
    char * input_tx_sam_file = argv[optind + 1];
    char * output_genome_sam_file = argv[optind + 2];

    FILE * input_tx_sam_fh = open_or_die(input_tx_sam_file, "r", "Input transcript-based SAM file");
    FILE * input_genome_sam_header_fh = open_if_present(input_genome_sam_header_file, "r");
    FILE * output_genome_sam_fh = open_or_die(output_genome_sam_file, "w", "Output genome-projected SAM file");


    // input SAM file is allowed to have or not have seq / qual.
    // basically, they are just payload, and unaffected by the transformation
    bool allow_absent_seq_qual = true;


    SamOrder genome_sam_order(SAM_POSITION_RID, "ALIGN");
    SamOrder tx_sam_order(SAM_POSITION_RID, "ALIGN");

    SAM_QNAME_FORMAT tx_qname_fmt = tx_sam_order.InitFromFile(input_tx_sam_fh);
    tx_sam_order.AddHeaderContigStats(input_tx_sam_fh);

    genome_sam_order.InitFromChoice(tx_qname_fmt);
    genome_sam_order.AddHeaderContigStats(input_genome_sam_header_fh);

    //here, construct a dummy low bound SamLine entry.
    
    SamLine::SetGlobalFlags(SAM_NON_INTERPRETED);
    char rname_buffer[256];
    std::fill(rname_buffer, rname_buffer + 256, 'A');

    SamLine output_lowbound(DATA_LINE, "output_bound", 0, rname_buffer, 0, 0, "", "*", 0, 0, "", "", "");

    SamLine::SetGlobalFlags(tx_qname_fmt);

    bool paired_reads_are_same_stranded = false;

    char const* species = "dummy";

    std::set<SequenceProjection> genome_to_tx = 
        gtf_to_sequence_projection(gtf_file, species);

    //needed if we are expanding reads that are originally mapped to a
    //transcript to genomic coordinates
    std::set<SequenceProjection> tx_to_genome;

    // std::set<cis::dna_t> contigs;

    // typedef std::set<cis::dna_t>::const_iterator DNAS_ITER;
    // DNAS_ITER dnas_iter;

    // for (SP_ITER gt_iter = genome_to_tx.begin();
    //      gt_iter != genome_to_tx.end(); ++gt_iter)
    // {
    //     contigs.insert(cis::dna_t((*gt_iter).species, (*gt_iter).source_dna,
    //                               0, INT64_MAX));
    // }


    //this maps tx_extent regions to sequence projections.  when in
    //expansion mode, the sequence projection will be
    //transcript->genome 't2g', and be iterators into tx_to_genome.
    //otherwise, these will be iterators into genome_to_tx
    
    PROJ_MAP tx_name_to_projection;
    

    std::map<cis::region const*, SP_ITER> tx_extent_to_projection;
    typedef std::map<cis::region const*, SP_ITER>::const_iterator MP_ITER;

    typedef std::pair<cis::REGIONS_MULTI::iterator, bool> REG_INS;

    // 1. tx_to_genome  (projections of transcripts to the genome)
    // 2. tx_name_to_projection.  (tx names pointing to tx_to_genome iterators)
    for (SP_ITER gt_iter = genome_to_tx.begin();
         gt_iter != genome_to_tx.end(); ++gt_iter)
    {
        SP_ITER tg_iter = tx_to_genome.insert(tx_to_genome.begin(), 
                                              InvertProjection((*gt_iter)));
        
        //transcript-to-genome projection for this transcript
        SequenceProjection const& tx_t2g_projection = *tg_iter;
        
        char * tx_name = new char[tx_t2g_projection.source_dna.size() + 1];
        strcpy(tx_name, tx_t2g_projection.source_dna.c_str());
        std::pair<NP_ITER, bool> ins = 
            tx_name_to_projection.insert(std::make_pair(tx_name, tg_iter));

        if (! ins.second)
        {
            SequenceProjection const& p1 = *(*ins.first).second;
            SequenceProjection const& p2 = (*tg_iter);

            fprintf(stderr, "Error: duplicate transcript named %s encountered.\n"
                    "First instance: %s\t%s\t%s\n"
                    "Second instance: %s\t%s\t%s\n", 
                    tx_name,
                    p1.source_dna.c_str(), (p1.same_strand ? "+" : "-"), Cigar::ToString(p1.cigar).c_str(),
                    p2.source_dna.c_str(), (p2.same_strand ? "+" : "-"), Cigar::ToString(p2.cigar).c_str());
            exit(1);
        }
    }

    //here need to determine / check consistency the genome contig
    //order implied by the transcriptome alignment
    bool consistent_order = 
        CheckProjectionOrder(genome_sam_order.contig_offsets,
                             tx_sam_order.contig_offsets,
                             tx_to_genome);

    if (! consistent_order)
    {
        fprintf(stderr, "Cannot continue.\n");
        exit(1);
    }


    // SP_ITER first_proj = tx_to_genome.end();
    // SP_ITER second_proj = tx_to_genome.end();

    char prev_rname[1024] = "";

    SamLine * samline;
    SamBuffer input_buffer(&tx_sam_order, paired_reads_are_same_stranded);
    SamBuffer output_buffer(&genome_sam_order, paired_reads_are_same_stranded);

    // dnas_iter = contigs.end();


    if (input_genome_sam_header_fh != NULL)
    {
        PrintSAMHeader(&input_genome_sam_header_fh, output_genome_sam_fh);
    }
    close_if_present(input_genome_sam_header_fh);

    char fake_samline[1024];

    size_t prev_pos_index = 0;
    size_t cur_pos_index = 0;

    std::pair<SamLine const*, bool> insert_result;

    CONTIG_OFFSETS::const_iterator contig_iter = input_buffer.sam_order->contig_offsets.begin();

    project_aux_data * thread_data_array = new project_aux_data[num_threads];

    fprintf(stderr, "Starting Projection...\n");
    fflush(stderr);

    for (size_t c = 0; c != num_threads; ++c)
    {
        project_aux_data * p = & thread_data_array[c];

        p->tx_name_to_projection = &tx_name_to_projection;
        p->inserts_are_introns = inserts_are_introns;
        p->add_cufflinks_xs_tag = add_cufflinks_xs_tag;
        p->record_ordering = output_buffer.sam_order;

    }


    while (! feof(input_tx_sam_fh))
    {
        // parse samline
        samline = new SamLine(input_tx_sam_fh, allow_absent_seq_qual);

        switch (samline->parse_flag)
        {
        case END_OF_FILE: 
            delete samline;
            ProjectTranscriptEntries(tx_name_to_projection, input_buffer,
                                     output_buffer, prev_rname,
                                     &output_lowbound,
                                     inserts_are_introns,
                                     add_cufflinks_xs_tag,
                                     num_threads,
                                     thread_data_array);
            break;

        case HEADER:
            delete samline;
            break;

        case PARSE_ERROR:
            fprintf(stderr, "Parse error in input sam file %s", input_tx_sam_file);
            exit(1);
            break;

        case DATA_LINE:

            //check well-orderedness of input
            sprintf(fake_samline, "fake\t%i\t%s\t%Zu", samline->flag, samline->rname, samline->pos);
            cur_pos_index = (input_buffer.sam_order->*(input_buffer.sam_order->sam_index))(fake_samline);
            if (cur_pos_index < prev_pos_index)
            {
                fprintf(stderr, "Error: input is not sorted by [rname, flag, pos] fields\n"
                        "Please sort by 'ALIGN' using align_eval sort\n");
                exit(1);
            }
            prev_pos_index = cur_pos_index;

            if (samline->query_unmapped())
            {
                samline->print(output_genome_sam_fh, true);
                delete samline;
                insert_result.second = false;
            }
            else
            {
                //initialize
                samline->SetFlattenedPosition(input_buffer.sam_order->contig_offsets, &contig_iter);
                insert_result = input_buffer.insert(samline);
            }

            if (insert_result.second && strcmp(prev_rname, samline->rname) != 0)
            {

                ProjectTranscriptEntries(tx_name_to_projection, input_buffer,
                                         output_buffer, prev_rname,
                                         &output_lowbound,
                                         inserts_are_introns,
                                         add_cufflinks_xs_tag,
                                         num_threads,
                                         thread_data_array);

                strcpy(prev_rname, samline->rname);
                output_buffer.purge(output_genome_sam_fh, NULL, NULL, &output_lowbound);
            } // if on new contig
            
        }
    }
    output_buffer.purge(output_genome_sam_fh, NULL, NULL, NULL);

    fclose(input_tx_sam_fh);
    fclose(output_genome_sam_fh);

    delete thread_data_array;

    for (NP_ITER tn = tx_name_to_projection.begin();
         tn != tx_name_to_projection.end(); ++tn)
    {
        delete (*tn).first;
    }
    
    for (MP_ITER tex = tx_extent_to_projection.begin();
         tex != tx_extent_to_projection.end(); ++tex)
    {
        delete (*tex).first;
    }

    return 0;
}
