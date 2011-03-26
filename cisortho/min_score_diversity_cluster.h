//defines the procedure for determining whether a certain
//neighborhood of hits constitutes a cluster based on having
//a minimum number of different hits at a minimum number of scores...
class min_score_diversity_cluster : public diversity_cluster {

	vector<vector<int> > minimum_scores;

 public:
	min_score_diversity_cluster(char const*, char const*, int);

	//bool trigger(RIT_M);
	void add_stat(RIT_M_R); //add the statistic of one hit to the cluster statistics
	void drop_stat(RIT_M_R); //drop the statistic of one hit

};

