#ifndef _ITEMGRAPH_H_
#define _ITEMGRAPH_H_

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <deque>
#include <utility>
#include "utils.h"
#include "sfmt/SFMT.h"

#define IF_TRACE(args) ;

namespace _Cide {
	
	class itemGraph {
		
	public: 
		
        //std::vector< std::vector< double> > probT;
        std::vector< std::vector<int> > hyperGT;
        
        // RR sets sampling related
        sfmt_t sfmtSeed{}; //struct
		std::vector<bool> visit; 
		std::vector<int> visit_mark;
		std::deque<int> q;
		
        void generateRRSample(vector<int> &targetNodes, int64 prevSize, int64 newSize);
        int BuildHypergraphNode(int uStart, int64 hyperiiid);
        
        itemGraph(int n); 
        ~itemGraph(void);
		
	};
	
}


#endif
