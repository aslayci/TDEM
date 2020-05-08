#ifndef _ITEMGRAPH_C_
#define _ITEMGRAPH_C_

#include "itemGraph.h"

namespace _Cide{
    
    itemGraph::itemGraph(int n) {
        visit_mark = std::vector<int>(n,0);
        visit = std::vector<bool>(n,false);
        
        sfmt_init_gen_rand(&sfmtSeed , 95082);
    }
    
    void itemGraph::generateRRSample(vector<int> &targetNodes, int64 prevSize, int64 newSize) {
        for (int64 z = prevSize; z < newSize; z++) {
            hyperGT.emplace_back(std::vector<int>());
            //hyperGT.push_back(std::vector<int>());
            BuildHypergraphNode(targetNodes[z], z);

        }
    }
    
    itemGraph::~itemGraph() {
        cout << "destructor called for itemGraph class" << endl;
    }

    //, std::vector< std::vector<int> > &hyperGT
    // only contain the node ids, rather than the pairs.
    int itemGraph::BuildHypergraphNode(int uStart, int64 hyperiiid) {
        int n_visit_edge=1;
        
        hyperGT[hyperiiid].push_back(uStart);
        
        int n_visit_mark = 0;
        
        q.clear();
        q.push_back(uStart);
        visit_mark[n_visit_mark++] = uStart; // repeatedly use this, not not initialize again.
        visit[uStart] = true;
        
        while(!q.empty()) {
            int expand = q.front();
            q.pop_front();
            
            int i = expand;
            for(int j = 0; j < (int)graphT[i].size(); j++){
                int v = graphT[i][j]; //parent of u in the original graph G
                n_visit_edge++;
                double randDouble = sfmt_genrand_real1(&sfmtSeed);
                if(randDouble > ((double)1.0/ (double)graphT[i].size())) // calculated on the fly
                    continue;
                if(visit[v])
                    continue;
                if(!visit[v]) {
                    visit_mark[n_visit_mark++] = v;
                    visit[v]=true;
                }
                q.push_back(v);
                hyperGT[hyperiiid].push_back(v);
            }
        }
        
        for(int i = 0; i < n_visit_mark; i++)
            visit[visit_mark[i]]=false;
        
        return n_visit_edge; 
    }
    
}

#endif

