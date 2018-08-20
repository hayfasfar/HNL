#ifndef HNL_HeavyNeutralLeptonAnalysis_BigNtuple
#define HNL_HeavyNeutralLeptonAnalysis_BigNtuple
/* 
	 Class: BigNtuple
	 Simple interface class to hide all the ROOT I/O from the plugin and make it more readable
*/

class TTree;
namespace edm {
	class EventID;
}

class BigNtuple {
public:
	BigNtuple(){} //default, empty constructor

	void set_evtinfo(TTree* tree);
	void fill_evtinfo(const edm::EventID& id);
	
private:
	unsigned int lumi_ = 0;
	unsigned int run_ = 0;
	unsigned long long evt_ = 0;
};

#endif
