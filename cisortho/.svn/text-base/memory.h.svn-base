#ifndef _MEMORY_H
#define _MEMORY_H

#include <set>

//allows storage and retrieval of data such that any attempt
//to store the same value again results in returning the reference
//to the stored object, rather than storing it a second time

//the thing needs to conform to the requirements of std::set, namely
//be copy-constructible and less-than comparable.

namespace memory {

	//problem here:  operator() makes the object a function of its parent object, and thus
	//a const qualifier must be used...
	//but that is to be expected, since it is in fact const...

template <typename T>
 class store {

	 static std::set<T> _store;
	 typename std::set<T>::iterator _ptr;

 public:
	 store<T>() : _ptr(_store.end()) { }

	 T & operator()() const {
		 if (_ptr == _store.end()) { throw "No stored item" ; }
		 else return const_cast<T&>(*_ptr);
	 }

	 T & operator()() {
		 if (_ptr == _store.end()) { throw "No stored item" ; }
		 else return const_cast<T&>(*_ptr);
	 }

	 T & operator()(T t) { 
		 _ptr = _store.insert(t).first; 
		 return const_cast<T&>(*_ptr); 
	 }
 };


 template <typename T> std::set<T> store<T>::_store;



template <typename T>
 class cstore {

	 static std::set<T> _cstore;
	 typename std::set<T>::const_iterator _ptr;

 public:
	 cstore<T>() : _ptr(_cstore.end()) { }

	 T const& operator()() const {
		 if (_ptr == _cstore.end()) { throw "No stored item" ; }
		 else return *_ptr;
	 }

	 T const& operator()(T const& t) { _ptr = _cstore.insert(t).first; return *_ptr; }
 };


 template <typename T> std::set<T> cstore<T>::_cstore;



} // end namespace memory


#endif // _MEMORY_H
