#include "observed_data.h"

void Observed_Data::init_data_index(ObservedData::DATA_COUNTS * data_counts)
{
    ObservedData::JPD::const_iterator iter;
    size_t jpd_index = 0;
    for (iter = (*data_counts).begin(); iter != (*data_counts).end(); ++iter)
    {
        ObservedData * mutable_datum =
            const_cast<ObservedData *>((*iter).first);
        mutable_datum->set_jpd_index(jpd_index);
        ++jpd_index;
    }
}


//checks one datum for whether it has nonzero probability in the jpd
size_t ObservedData::valid_counts_data(ObservedData::JPD const& prior_jpd_data)
{

    size_t number_invalid = 0;
    ObservedData::JPD::const_iterator iter;
    for (iter = this->raw_counts_data.begin();
         iter != this->raw_counts_data.end();
         ++iter)
    {
        double jpd_slice_prob = 0.0;
        for (size_t b = 0; b != NBASES; ++b)
        {
            jpd_slice_prob += prior_jpd_data[b][datum->jpd_index()];
        }
        number_invalid += jpd_slice_prob > 0.0 ? 1 : 0;
    }
    return number_invalid;
}
