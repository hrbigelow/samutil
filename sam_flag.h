#ifndef _SAM_FLAG_H
#define _SAM_FLAG_H

#include <cstddef>

struct SamFlag
{
    unsigned int multi_fragment_template : 1;     // SAM
    unsigned int all_fragments_mapped : 1;        // SAM and rSAM
    unsigned int this_fragment_unmapped : 1;      // SAM
    unsigned int next_fragment_unmapped : 1;      // SAM
    unsigned int this_fragment_on_neg_strand : 1; // SAM
    unsigned int next_fragment_on_neg_strand : 1; // SAM
    unsigned int first_fragment_in_template : 1;  // SAM
    unsigned int last_fragment_in_template : 1;   // SAM
    unsigned int alignment_not_primary : 1;       // SAM and rSAM
    unsigned int failed_quality_check : 1;        // SAM and rSAM
    unsigned int pcr_or_optical_duplicate : 1;    // SAM and rSAM
    unsigned int template_layout : 1;             //         rSAM
    unsigned int : 4; // padding to 16 bits

    unsigned int is_rsam_format : 8;              //         rSAM
    unsigned int num_fragments_in_template : 8;   //         rSAM
    unsigned int read_layout : 32;                //         rSAM

    size_t get_raw() const;
    void set_raw(size_t raw);

};



#endif // _SAM_FLAG_H
