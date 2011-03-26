#include <vector>
#include <iostream>
#include "boost/function.hpp"
//#include "/razor/5/users/bigelow/src/FC++.1.5/prelude.h"


using std::cerr;
using std::endl;
using std::cout;
using boost::function;

template <typename Element, typename Seed>
Seed Fold(Seed seed, function<Seed (Seed &, Element *)> fold_fun, function<Element * (void)> next){
	Element *cur = next();
	if (cur == NULL){ return seed; }
	else {
		Seed newseed = fold_fun(seed, cur);
		return Fold(newseed, fold_fun, next);
	}
}


/* template <typename Element, typename Seed> */
/* Seed Fold2(Seed seed, fcpp::Fun2<Seed, Element *, Seed> fold_fcn, fcpp::Fun0<Element *> next){ */
/* 	Element *cur = next(); */
/* 	if (cur == NULL){ return seed; } */
/* 	else { */
/* 		Seed newseed = fold_fcn(seed, cur); */
/* 		return Fold2(newseed, fold_fcn, next); */
/* 	} */
/* } */


typedef std::vector<int> coord;

enum Exceptions {
	NoMoreElements
};




void increment_counter_aux(const coord & sizes, coord & ctr, int d) throw (Exceptions) {
	
	if (d == (int)ctr.size()) throw NoMoreElements;
	
	ctr[d]++;
	if (ctr[d] == sizes[d]){
		ctr[d] = 0;
		increment_counter_aux(sizes, ctr, d+1);
	}
}


void increment_counter(const coord & sizes, coord & counter){
	increment_counter_aux(sizes, counter, 0);
}


//how should we handle zero-sized dimensions?
//do not handle them...
template <typename Data> struct CartesianProductEnum {
	
	vector<int> counter;
	vector<int> sizes;
	vector<Data> current;
	vector<vector<Data> >source;
	
	CartesianProductEnum(vector<vector<Data> >& matrix, vector<Data> &initial)  {
		
		current = initial;
		source = matrix;
		if (current.size() != source.size()) current.resize(source.size());
		
		counter.resize(matrix.size());
		sizes.resize(matrix.size());
		fill(counter.begin(), counter.end(), 0);
		if (! counter.empty()) counter[0] = -1;
		for (unsigned int i=0; i < source.size(); i++) {
			sizes[i] = source[i].size();
			if (sizes[i] == 0) {
				cerr<<"CartesianProductEnum: Cannot handle zero-sized dimensions."<<endl;
				exit(75);
			}
		}
/* 		cout<<"sizes:"; */
/* 		for (unsigned int i=0; i < counter.size(); i++) cout<<'\t'<<sizes[i]; */
/* 		cout<<endl; */
	}
	
	vector<Data> * Next(){

		try { increment_counter(sizes, counter); }
		catch (Exceptions e){ 
/* 			cout<<"END-------"<<endl<<endl; */
			return NULL; 
		}

/* 		cout<<"count:"; */
/* 		for (unsigned int i=0; i < counter.size(); i++) cout<<'\t'<<counter[i]; */
/* 		cout<<endl; */


		for(unsigned int i = 0; i < counter.size(); i++)
			current[i] = source[i][counter[i]];

		return &current;
	}
};
		
