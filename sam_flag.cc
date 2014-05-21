#include "sam_flag.h"


size_t SamFlag::get_raw() const
{
    size_t raw = static_cast<size_t>
        (this->multi_fragment_template % 2
         | (this->all_fragments_mapped % 2)<<1
         | (this->this_fragment_unmapped % 2)<<2
         | (this->next_fragment_unmapped % 2)<<3
         | (this->this_fragment_on_neg_strand % 2)<<4
         | (this->next_fragment_on_neg_strand % 2)<<5
         | (this->first_fragment_in_template % 2)<<6
         | (this->last_fragment_in_template % 2)<<7
         | (this->alignment_not_primary % 2)<<8
         | (this->failed_quality_check % 2)<<9
         | (this->pcr_or_optical_duplicate % 2)<<10
         | (this->template_layout % 2)<<11
         )
        | static_cast<size_t>(this->is_rsam_format % (1<<8))<<16
        | static_cast<size_t>(this->num_fragments_in_template % (1<<8))<<24
        | static_cast<size_t>(this->read_layout % (1L<<32))<<32;

    return raw;
}


void SamFlag::set_raw(size_t raw)
{
    this->multi_fragment_template = static_cast<unsigned int>(raw % 2);
    this->all_fragments_mapped = static_cast<unsigned int>((raw>>1) % 2);
    this->this_fragment_unmapped = static_cast<unsigned int>((raw>>2) % 2);
    this->next_fragment_unmapped = static_cast<unsigned int>((raw>>3) % 2);
    this->this_fragment_on_neg_strand = static_cast<unsigned int>((raw>>4) % 2);
    this->next_fragment_on_neg_strand = static_cast<unsigned int>((raw>>5) % 2);
    this->first_fragment_in_template = static_cast<unsigned int>((raw>>6) % 2);
    this->last_fragment_in_template = static_cast<unsigned int>((raw>>7) % 2);
    this->alignment_not_primary = static_cast<unsigned int>((raw>>8) % 2);
    this->failed_quality_check = static_cast<unsigned int>((raw>>9) % 2);
    this->pcr_or_optical_duplicate = static_cast<unsigned int>((raw>>10) % 2);
    this->template_layout = static_cast<unsigned int>((raw>>11) % 2);

    this->is_rsam_format = static_cast<unsigned int>((raw>>16) % (1<<8));
    this->num_fragments_in_template = static_cast<unsigned int>((raw>>24) % (1<<8));
    this->read_layout = static_cast<unsigned int>((raw>>32) % (1<<8));

}
