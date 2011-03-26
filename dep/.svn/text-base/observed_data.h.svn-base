class ObservedDatum {

    size_t _jpd_index; //position in the JPD array this observed data will be placed
    size_t jpd_index() const { return _jpd_index; }
    void set_jpd_index(size_t i) { _jpd_index = i; }
    virtual ObservedData() = 0;
    virtual ~ObservedData();
};



class ObservedData {

 public:

};



//base class for overloading the process of parsing the observed data
//and expanding it from raw form to jpd with the founder base identity
//provides the raw counts, the jpd, and the base marginal
class ObservedData {

 public:
    typedef double * JPD[NBASES];
    typedef std::map<ObservedDatum *, size_t> DATA_COUNTS;
    
    JPD complete_jpd;
    
    double founder_base_marginal[NBASES];
    
    DATA_COUNTS raw_counts_data;

    //loads prior data from a tabular file of counts, expands it into
    //a JPD including the founder base identity as an extra dimension
    virtual void load_from_summary_file(char const* data_file);

    virtual void expand_to_jpd();

    virtual void load_locus(void * locus_data, void * extra_params);

    void init_jpd_index(JPD * jpd_data);
    
    size_t valid_counts_data(JPD const& prior_data);


    virtual ObservedData() = 0;
    virtual ~ObservedData() = 0; //clean up
};
