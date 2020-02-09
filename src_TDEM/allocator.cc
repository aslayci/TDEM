#include "allocator.h"
#include "anyoption.h"
#include "iheap.h"
#include <ctime>
#include <numeric>

namespace _Cide{
    
    allocator::allocator(AnyOption* opt1) {
        
        opt = opt1;
        delim = " \t";
        
        sfmt_init_gen_rand(&sfmtSeed , 95082);
        
        n = strToInt(opt->getValue("n"));
        m = strToInt(opt->getValue("m"));
        k = strToInt(opt->getValue("k"));
        k_r = strToInt(opt->getValue("kr"));
        k_b = strToInt(opt->getValue("kb"));
        tao = (int) k_b/k_r;
        P = Permutation(n, n-k_r*(tao + 1)) / (double) (pow(fact(tao), k_r) * fact(k_r));

        // TO ensure that P is correct
        cout << "tao is " << tao << endl;
        cout << "P is " << P << endl;
        cout << "Compare is " << pow(2.71, logcnk(nrPairs, k)) << endl;

        nrItems = strToInt(opt->getValue("nrItems"));
        epsilon = strToDouble(opt->getValue("epsilon"));
        ell = strToDouble(opt->getValue("ell"));
        outFolderName = opt->getValue("outputFolder");
    
        nrPairs = n * nrItems; // multiple items
        
        for(int i = 0; i < n; i++)
            graphT.emplace_back(vector<int>());

        prevSize = 0; // when sample the first node
        for(int i = 0; i < nrPairs; i++)
            hyperG.emplace_back(std::vector<int>());

        // create new pairs-related functions
        hyper_degreepairs = std::vector<double>(n*n, 0);
        for (int i = 0; i < n*n; ++i) {
           hyperGpairs.emplace_back(std::vector<int>());
        }

        rcList = new itemGraphList();
        //item-specific probT is aligned with graphT for each item
        for(int i = 0; i < nrItems; i++) {
            //_Cide::itemGraph *ig = new _Cide::itemGraph(n);
            auto *ig = new _Cide::itemGraph(n);
            for (int j = 0; j < n; j++)
                //ig->probT.push_back(std::vector< double>());
                ig->probT.emplace_back(std::vector< double>());
            rcList->push_back(ig);
        }

        nodeDegree = std::vector<double>(n,0.0);
        
        cout << "nr items " << nrItems << endl;
        cout << "assignment size " << k << endl;
        
        readNodeLeaningsFile();
        readItemLeaningsFile();
        
        readTICGraph();
        cout << "successfully read everything " << endl;
        arrangeOutputFiles();
        cout << "starting assignment for tdem" << endl;
        tdem();
        
// comment out below to run the baselines
//        degreeClose(theta, k);
//        degreeFar(theta, k);
//        degreeWeighted(theta, k);
        
    }
    
    double allocator::tdem() {
        
        cout << "computation for lower bounding OPT started " << endl;

        clock_t common_begin = clock();
        double lb = lowerBoundOPT();
        cout << "lower bound identified: " << lb << endl;
        
        theta = (8 + 4/3.0 * epsilon) * n * (ell*log(n) + log(1 + log2(n)) +  log(P)) / (epsilon * epsilon * lb);
        
        cout << "final sample size " << theta << endl;
        generateRCSets(theta);
        
        clock_t common_end = clock();
        duration_common = double(common_end - common_begin) / CLOCKS_PER_SEC;
        
        double greedySolution = n * rcGreedy(max(prevSize, theta), k, true);
        return greedySolution; 
    }

    
    double allocator::lowerBoundOPT() {
        
        double epsilon_1 = epsilon;
        
        for (int i = 1; i < log2(n); i++) {
//            cout << "here x = " << x << endl;
            float x = n / (double) pow(2.0, i);
            int64 theta_x = (2+2/3.0 * epsilon_1)* (ell*log(n) + log(P) + log(log2(n) + 1)) *n / (x * epsilon_1 * epsilon_1);
            cout << "here theta_x " << theta_x << endl;
            generateRCSets(theta_x);

            double ept = rcGreedy(theta_x, k, false);
            cout << "ept " << ept << endl;

            if (ept*n > ((1+epsilon_1) * x)) {
                double lowerBound = ept * n / (1 + epsilon_1);
                return lowerBound;
            }
        }
        cout << "returning naive lower bound  " << endl;
        double naive = (double) k * 2.0;
        return naive;
        
    }

    void allocator::generateRCSets(int64 newSize) {

        if(newSize < prevSize){
            return;
        }
        // sample target nodes
        for (int i = prevSize; i < newSize; i++) {
            int randTarget = sfmt_genrand_uint32(&sfmtSeed) % n;
            targetNodes.push_back(randTarget);
        }
        
        // expand coordinated RR sets samples of items
        for(int itemID = 0; itemID < nrItems; itemID++) {
            rcList->at(itemID)->generateRRSample(targetNodes, prevSize, newSize);
        }

        // rcList -> itemID -> hyperGT[rcID]: store a list of the nodes.
        double tmpsize = createhyperGTPairs(newSize);
        cout << "size of hyPerGTPairs is " << tmpsize << endl;

        for (int rcID = prevSize; rcID < newSize; ++rcID) {
            int tpsize = hyperGTpairs[rcID].size();
            for (int i = 0; i < tpsize; ++i) {
                int x = hyperGTpairs[rcID][i]; // because x is too large.
                hyperGpairs[x].push_back(rcID);
                hyper_degreepairs[x] += 1;
            }
            int tprsize = rcList->at(0)->hyperGT[rcID].size();
            for (int j = 0; j < tprsize; ++j) {
                int y = rcList->at(0)->hyperGT[rcID][j];
                hyperG[y].push_back(rcID); // careful about the repeatation
            }
            int tpbsize = rcList->at(1)->hyperGT[rcID].size();
            for (int j = 0; j < tpbsize; ++j) {
                int y = rcList->at(1)->hyperGT[rcID][j];
                hyperG[y+n].push_back(rcID); // careful about the repeatation
            }
        }

        prevSize = newSize;

    }

    double allocator::createhyperGTPairs(int64 newSize) {
        // use to create reverse rr pairs
        //update hyperG and hyper_degree
        cout << "newSize is " << newSize << endl;
        for (int z = prevSize; z < newSize; ++z) {
            hyperGTpairs.emplace_back(std::vector<int>()); // init
        }

        for (int rcID = prevSize; rcID < newSize; ++rcID) {
            for (int i = 0; i < rcList->at(0)->hyperGT[rcID].size(); ++i) {
                for (int j = 0; j < rcList->at(1)->hyperGT[rcID].size(); ++j) {
                   int t = rcList->at(0)->hyperGT[rcID][i]*n + rcList->at(1)->hyperGT[rcID][j];
                   hyperGTpairs[rcID].push_back(t);
                }
            }
        }

        return hyperGTpairs.size();
    }


    double allocator::rcGreedy(int64 rcSampleSize, int k, bool extraResults) {

        // For initialization the k_r and k_b, to delete later.

        clock_t begin = clock();
        std::vector<double> ExpScore(n*n, 0); // use to store hyper_degreepairs

        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (i != j){
                    ExpScore[i*n + j] = hyper_degreepairs[i*n + j];
                }
            }
        }

        for (int k1 = 0; k1 < n; ++k1) {
            ExpScore[k1*n + k1]  = 0;
        }

        double totalExpScore = 0;
        int bestPairID, bestRedNodeID, bestBlueNodeID, bestItemID;
        std::vector<int> isValidRed = std::vector<int>(n,tao);
        std::vector<int> isValidBlue = std::vector<int>(n,1);
        std::vector<bool> isCovered = std::vector<bool>(rcSampleSize, false);


        seedSet.clear();
        seedScores.clear(); // seedScores should be cleared as well.
        while ((int)seedSet.size() < k_b) {
            // here I use the vector to store the results, but can use the priority queue.
            bestPairID = std::max_element(ExpScore.begin(), ExpScore.end()) - ExpScore.begin();
            //cout << "ExpScore is " << ExpScore[bestPairID] << "rc sample size is " <<rcSampleSize << " n is " << n << endl;
            bestRedNodeID = (int) bestPairID / n;
            bestBlueNodeID = (int) bestPairID % n;

            // if can not select, then just pop out
            if (isValidRed[bestRedNodeID] == 0 || isValidBlue[bestBlueNodeID] == 0) {
                ExpScore[bestPairID] = 0;
                continue;
            }

            // push
            seedSet.push_back(bestPairID); // store
            //cout << "ExpScore is " << ExpScore[bestPairID] << "rc sample size is " <<rcSampleSize << "n is " << n << endl;
            seedScores.push_back(ExpScore[bestPairID] / (double) rcSampleSize); // store score

            // make update
            isValidRed[bestRedNodeID] -= 1;
            isValidRed[bestBlueNodeID] = 0;
            isValidBlue[bestBlueNodeID] = 0;
            isValidBlue[bestRedNodeID] = 0;

            if(isValidRed[bestRedNodeID] <= 0){
                for (int i = 0; i < n; ++i) {
                   ExpScore[bestRedNodeID*n + i] = 0;
                }
            }
            for (int i = 0; i < n; ++i) {
               ExpScore[i*n + bestBlueNodeID]  = 0;
               ExpScore[i*n + bestRedNodeID]  = 0;
               ExpScore[bestBlueNodeID*n + i] = 0;
            }
            

            // update the ExpScore of the pairs whose marginal gain change due to the best pair selection
            int tmpsize = hyperGpairs[bestPairID].size();
            for (int j = 0; j < tmpsize; j++) {
                int rcID = hyperGpairs[bestPairID][j]; // get the rcID
                if(isCovered[rcID] == false){ // ensure not to repeatedly --
                    int ttmpsize = hyperGTpairs[rcID].size();
                    for (int z = 0; z < ttmpsize; z++) {
                        int pairID = hyperGTpairs[rcID][z];
                        if (ExpScore[pairID] != 0) {
                            ExpScore[pairID]--;
                        }
                    }
                    isCovered[rcID] = true;
                }
           }
        } // end of while k loop

        for (int i1 = 0; i1 < seedScores.size(); ++i1) {
           totalExpScore += seedScores[i1];
        }

        ///////////// produce the additional results here
        ////////////  Remember should write the f function results not the g function result.
        if(extraResults) {
            
            clock_t end = clock();
            totalDuration = (double(end - begin) / CLOCKS_PER_SEC) + duration_common; // in seconds
            totalMemory = getCurrentMemoryUsage(); // in MB
            
            cout << "total tdem time taken (in seconds) " << totalDuration << " total memory (in mb) " << totalMemory << endl;
            
            double totMGs = 0;
            for (int i = 0; i < seedSet.size(); i++) {
                int nodeRed = (int) seedSet[i] / n;
                int nodeBlue = (int) seedSet[i] % n;
                totMGs += seedScores[i];
                writeInMasterOutputFile(nodeRed,nodeBlue,seedScores[i],totMGs,totalDuration,totalMemory);
            }
            cout << "total tdem div. exp. score " << totMGs << endl;
            
            // compute for each seed pair its prob. of reaching target nodes
            
            string extraFName1 = outFolderName + OS_SEP + "pairReachProbs_tdem.txt";
            ofstream extraStream1;
            
            if(extraStream1.is_open())
                extraStream1.close();
            
            extraStream1.open(extraFName1.c_str());
            
            if (extraStream1.is_open() == false) {
                cout << "Can't open file " << extraFName1  << " for writing" << endl;
                exit(1);
            }

            extraStream1 << "RednodeID" << " " << "BluenodeID" << " " <<  "targetNodeID" << " " << "reachProb" << endl;
            
            extraStream1 << std::fixed;
            extraStream1.precision(8);
            
            std::vector< vector<double> > pairNodeReachProbs;

            for (int i = 0; i < k_b; i++) {
                pairNodeReachProbs.push_back(std::vector<double>(n,0.0));
            }
            
            for (int i = 0; i < k_b; i++) {
                int seedPairID = seedSet[i];
                vector<int> rcsCovered = hyperGpairs[seedPairID]; // ids of RC sets that contain this best pair
                for (int j = 0; j < rcsCovered.size(); j++) {
                    int rcID = rcsCovered[j];
                    int target = targetNodes[rcID];
                    double tempVal = pairNodeReachProbs[i][target];
                    pairNodeReachProbs[i][target] = tempVal + 1.0;
                }
            }
            
            for (int i = 0; i < k_b; i++) {
                int seedPairID = seedSet[i];
                int RednodeID = (int) seedPairID / n;
                int BluenodeID = (int) seedPairID % n;
                for (int j = 0; j < n; j++) {
                    double tempVal = pairNodeReachProbs[i][j];
                    pairNodeReachProbs[i][j] = (double) n * tempVal / (double) rcSampleSize;
                    extraStream1 << RednodeID << " " << BluenodeID << " " <<  j << " " << pairNodeReachProbs[i][j] << endl;
                }
            }
            
            extraStream1.close();
            
            // compute f(seedset) for each target node v, should consider the other pairs.
            string extraFName2 = outFolderName + OS_SEP + "nodeExpLevels_tdem.txt";
            ofstream extraStream2;
            
            if(extraStream2.is_open())
                extraStream2.close();
            
            extraStream2.open(extraFName2.c_str());

            if (extraStream2.is_open() == false) {
                cout << "Can't open file " << extraFName2  << " for writing" << endl;
                exit(1);
            }
            
            extraStream2 << "targetNodeID" << " " << "final_score_f" << " " << "in total"<< endl;

            vector<double> nodeSampleExpLevels(rcSampleSize,0.0);
            vector<double> nodeExpLevels(n,0.0);
            std::vector<int> rrRed;
            std::vector<int> rrBlue;

            for (int i = 0; i < k_b; ++i) {
                int SeedPairID = seedSet[i];
                for (int j = 0; j < hyperG[SeedPairID / n].size(); ++j) {
                    int id = hyperG[SeedPairID / n][j];
                    rrRed.push_back(id);
                }
                for (int l = 0; l < hyperG[SeedPairID % n + n].size(); ++l) {
                    int id = hyperG[SeedPairID % n + n][l];
                    rrBlue.push_back(id);
                }
            }

            std::sort(rrRed.begin(), rrRed.end());
            std::sort(rrBlue.begin(), rrBlue.end());

            std::vector<int> rrRedBlue(rrRed.size() + rrBlue.size(), 0);
            set_intersection(rrRed.begin(), rrRed.end(), rrBlue.begin(), rrBlue.end(), rrRedBlue.begin());

            for (int zz = 0; zz < rrRedBlue.size(); ++zz) {
                int tmp = targetNodes[rrRedBlue[zz]];
                nodeExpLevels[tmp] += 1;
            }

            float totalscore = 0;
            
//            double kontrolAw = 0;
            for (int i = 0 ; i < n ; i++) {
                double tempVal = nodeExpLevels[i];
                nodeExpLevels[i] = (double) tempVal / (double) rcSampleSize;
                totalscore += nodeExpLevels[i];
//                kontrolAw += nodeExpLevels[i];
                extraStream2 << i << " " << nodeExpLevels[i] << " " << totalscore << endl;
            }

            extraStream2.close();
//
            cout << "Total score is " << totalscore << endl;
            
        } // end of extra results
        
        return totalExpScore;
    }
    
    double allocator::degreeClose(int64 rcSampleSize, int k) {
        
        clock_t begin = clock();
        
        priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond>heap;
        std::vector<double> divExpScore(nrPairs, 0);

        std::vector<double> rcMinSeen(rcSampleSize,-5.0);
        std::vector<double> rcMaxSeen(rcSampleSize, 5.0);
        
        vector<bool> isCovered(rcSampleSize, false);
        vector<bool> pairMark(nrPairs, true);
        
        double totalExpScore = 0;
        int bestPairID, bestNodeID, bestItemID;
        
        divExpScore = hyper_degree;
        
        for (int i = 0; i < n; i++) {
            double val = (double) nodeDegree[i];
            std::pair<int, double>pairVal(std::make_pair(i, val));
            heap.push(pairVal);
        }
        seedSet.clear();
        seedScores.clear();
        
        while ((int)seedSet.size() < k) {
            
            pair<int, double> pairVal = heap.top();
            bestNodeID = pairVal.first;

            double itMax = 2.0;
            for (int itemID = 0; itemID < nrItems; itemID++) {
                int pairIDtemp = bestNodeID + n * itemID;
                if (pairMark[pairIDtemp]) { // if this item is not assigned to this node yet
                    if(itMax > std::abs(nodeLeanings[bestNodeID] - itemLeanings[itemID])) {
                        itMax = std::abs(nodeLeanings[bestNodeID] - itemLeanings[itemID]);
                        bestItemID = itemID;
                    }
                }
            }
            
            bestPairID = bestNodeID + n * bestItemID;
            totalExpScore += divExpScore[bestPairID];
            seedSet.push_back(bestPairID);
            seedScores.push_back(n * divExpScore[bestPairID] / (double) rcSampleSize);
            pairMark[bestPairID] = false;
            
            // update the hyper_degree of the pairs whose marginal gain change due to the best pair selection
            vector<int> rcsAffected = hyperG[bestPairID]; // ids of RC sets that contain this best pair
            for (int j = 0; j < rcsAffected.size(); j++) {
                int rcID = rcsAffected[j];
                int target = targetNodes[rcID];
                double prevGain, currGain;
                
                // if it is the first time this RC is covered by a seed pair -- below the marginal gain update operations for this case
                if (!isCovered[rcID]) {
                    isCovered[rcID] = true;
                    rcMinSeen[rcID] = std::min(nodeLeanings[target], itemLeanings[bestItemID]);
                    rcMaxSeen[rcID] = std::max(nodeLeanings[target], itemLeanings[bestItemID]);
                    
                    // compute the amount of change in marginal gain for each item contain in this RC set
                    for(int itemID = 0; itemID < nrItems; itemID++) {
                        
                        prevGain = std::abs(nodeLeanings[target] - itemLeanings[itemID]); // mg to empty set
                        if (itemLeanings[itemID] >= rcMinSeen[rcID] &&  itemLeanings[itemID] <= rcMaxSeen[rcID]) {
                            currGain = 0;
                        }
                        else if(itemLeanings[itemID] < rcMinSeen[rcID]) {
                            currGain = rcMinSeen[rcID] - itemLeanings[itemID];
                        }
                        else if(itemLeanings[itemID] > rcMaxSeen[rcID]) {
                            currGain = itemLeanings[itemID] - rcMaxSeen[rcID];
                        }
                        vector<int> rrTemp = rcList->at(itemID)->hyperGT[rcID];
                        for (int z = 0; z < rrTemp.size(); z++) {
                            int pairID = rrTemp[z] + n * itemID;
                            if(pairMark[pairID]) {
                                double tempD = divExpScore[pairID];
                                divExpScore[pairID] = tempD - prevGain + currGain;
                                if(divExpScore[pairID] < 0) {
                                    divExpScore[pairID] = 0;
                                }
                            }
                        }
                    }
                }
                
                // if it was already covered before
                else {
                    if(itemLeanings[bestItemID] >= rcMinSeen[rcID] && itemLeanings[bestItemID] <= rcMaxSeen[rcID])
                        continue;
                    
                    double prevMin = rcMinSeen[rcID];
                    double prevMax = rcMaxSeen[rcID];
                    rcMinSeen[rcID] = std::min(prevMin, itemLeanings[bestItemID]);
                    rcMaxSeen[rcID] = std::max(prevMax, itemLeanings[bestItemID]);
                    
                    for(int itemID = 0; itemID < nrItems; itemID++) {
                        if (itemLeanings[itemID] >= prevMin &&  itemLeanings[itemID] <= prevMax) {
                            prevGain = 0;
                        }
                        else if(itemLeanings[itemID] < prevMin) {
                            prevGain = prevMin - itemLeanings[itemID];
                        }
                        else if(itemLeanings[itemID] > prevMax) {
                            prevGain = itemLeanings[itemID] - prevMax;
                        }
                        
                        if (itemLeanings[itemID] >= rcMinSeen[rcID] &&  itemLeanings[itemID] <= rcMaxSeen[rcID]) {
                            currGain = 0;
                        }
                        else if(itemLeanings[itemID] < rcMinSeen[rcID]) {
                            currGain = rcMinSeen[rcID] - itemLeanings[itemID];
                        }
                        else if(itemLeanings[itemID] > rcMaxSeen[rcID]) {
                            currGain = itemLeanings[itemID] - rcMaxSeen[rcID];
                        }
                        
                        vector<int> rrTemp = rcList->at(itemID)->hyperGT[rcID];
                        for (int z = 0; z < rrTemp.size(); z++) {
                            int pairID = rrTemp[z] + n * itemID;
                            if(pairMark[pairID]) {
                                double tempD = divExpScore[pairID];
                                divExpScore[pairID] = tempD - prevGain + currGain;
                                if(divExpScore[pairID] < 0) {
                                    divExpScore[pairID] = 0;
                                }
                            }
                        }
                    }
                } // end-if this RC was already covered before
                
            }// end of for loop for affected
        } // end of while k loop
        
        
        
        clock_t end = clock();
        totalDuration = (double(end - begin) / CLOCKS_PER_SEC) + duration_common; // in seconds
        totalMemory = getCurrentMemoryUsage(); // in MB
        
        cout << "total degree-close time taken (in seconds) " << totalDuration << " total memory (in mb) " << totalMemory << endl;
        
        string extraFName0 = outFolderName + OS_SEP + "assignment_degreeClose.txt";
        ofstream extraStream0;
        if(extraStream0.is_open())
            extraStream0.close();
        
        extraStream0.open(extraFName0.c_str());
        
        if (extraStream0.is_open() == false) {
            cout << "Can't open file " << extraFName0  << " for writing" << endl;
            exit(1);
        }
        
        extraStream0 << "nodeID" << " " << "itemID" << " " <<  "mgScore" << " " << "cumScore" << " " << "duration(sec)" << " " << "memory(mb)" << endl;
        double totMGs = 0;
        for (int i = 0; i < seedSet.size(); i++) {
            int node = (int) seedSet[i] % n;
            int item = (int) seedSet[i] / n;
            totMGs += seedScores[i];
            extraStream0 << node << " " << item << " " << seedScores[i] << " " << totMGs << " " << totalDuration << " " << totalMemory << endl;
        }
        cout << "total degree-close div. exp. score " << totMGs << endl;
        
        extraStream0.close();
        
        string extraFName1 = outFolderName + OS_SEP + "pairReachProbs_degreeClose.txt";
        ofstream extraStream1;
        
        if(extraStream1.is_open())
            extraStream1.close();
        
        extraStream1.open(extraFName1.c_str());
        
        if (extraStream1.is_open() == false) {
            cout << "Can't open file " << extraFName1  << " for writing" << endl;
            exit(1);
        }
        
        extraStream1 << "nodeID" << " " << "itemID" << " " <<  "targetNodeID" << " " << "reachProb" << endl;
        
        extraStream1 << std::fixed;
        extraStream1.precision(8);
        
        vector< vector<double> > pairNodeReachProbs;
        
        for (int i = 0; i < k; i++) {
            pairNodeReachProbs.push_back(std::vector<double>(n,0.0));
        }
        
        for (int i = 0; i < k; i++) {
            int seedPairID = seedSet[i];
            vector<int> rcsCovered = hyperG[seedPairID]; // ids of RC sets that contain this best pair
            for (int j = 0; j < rcsCovered.size(); j++) {
                int rcID = rcsCovered[j];
                int target = targetNodes[rcID];
                double tempVal = pairNodeReachProbs[i][target];
                pairNodeReachProbs[i][target] = tempVal + 1.0;
            }
        }
        
        for (int i = 0; i < k; i++) {
            int seedPairID = seedSet[i];
            int nodeID = (int) seedPairID % n;
            int itemID = (int) seedPairID / n;
            for (int j = 0; j < n; j++) {
                double tempVal = pairNodeReachProbs[i][j];
                pairNodeReachProbs[i][j] = (double) n * tempVal / (double) rcSampleSize;
                extraStream1 << nodeID << " " << itemID << " " <<  j << " " << pairNodeReachProbs[i][j] << endl;
            }
        }
        
        extraStream1.close();
        
        // compute fv(seedset) for each target node v
        string extraFName2 = outFolderName + OS_SEP + "nodeExpLevels_degreeClose.txt";
        ofstream extraStream2;
        
        if(extraStream2.is_open())
            extraStream2.close();
        
        extraStream2.open(extraFName2.c_str());
        
        if (extraStream2.is_open() == false) {
            cout << "Can't open file " << extraFName2  << " for writing" << endl;
            exit(1);
        }
        
        extraStream2 << "targetNodeID" << " " << "exposure_level" << endl;
        vector<double> nodeExpLevels(n,0.0);
        
        for (int rcID = 0; rcID < rcSampleSize; rcID++) {
            int target = targetNodes[rcID];
            if(isCovered[rcID]) {
                nodeExpLevels[target] += (rcMaxSeen[rcID] - rcMinSeen[rcID]);
            }
        }
        
        for (int i = 0 ; i < n ; i++) {
            double tempVal = nodeExpLevels[i];
            nodeExpLevels[i] = (double) n * tempVal / (double) rcSampleSize;
            extraStream2 << i << " " << nodeExpLevels[i] << endl;
        }
        
        extraStream2.close();
        
        return totalExpScore / (double) rcSampleSize;
    }

    double allocator::degreeFar(int64 rcSampleSize, int k) {
        
        clock_t begin = clock();
        
        priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond>heap;
        std::vector<double> divExpScore(nrPairs, 0);

        std::vector<double> rcMinSeen(rcSampleSize,-5.0);
        std::vector<double> rcMaxSeen(rcSampleSize, 5.0);
        
        vector<bool> isCovered(rcSampleSize, false);
        vector<bool> pairMark(nrPairs, true);
        
        double totalExpScore = 0;
        int bestPairID, bestNodeID, bestItemID;
        
        
        divExpScore = hyper_degree;
        
        for (int i = 0; i < n; i++) {
            double val = (double) nodeDegree[i];
            std::pair<int, double>pairVal(std::make_pair(i, val));
            heap.push(pairVal);
        }
        
        seedSet.clear();
        seedScores.clear();
        
        while ((int)seedSet.size() < k) {
            
            pair<int, double> pairVal = heap.top();
            bestNodeID = pairVal.first;

            double itMax = 0;
            for (int itemID = 0; itemID < nrItems; itemID++) {
                int pairIDtemp = bestNodeID + n * itemID;
                if (pairMark[pairIDtemp]) { // if this item is not assigned to this node yet
                    if(itMax < std::abs(nodeLeanings[bestNodeID] - itemLeanings[itemID])) {
                        itMax = std::abs(nodeLeanings[bestNodeID] - itemLeanings[itemID]);
                        bestItemID = itemID;
                    }
                }
            }
            
            bestPairID = bestNodeID + n * bestItemID;
            totalExpScore += divExpScore[bestPairID];
            seedSet.push_back(bestPairID);
            seedScores.push_back(n * divExpScore[bestPairID] / (double) rcSampleSize);
            pairMark[bestPairID] = false;
            
            // update the hyper_degree of the pairs whose marginal gain change due to the best pair selection
            vector<int> rcsAffected = hyperG[bestPairID]; // ids of RC sets that contain this best pair
            for (int j = 0; j < rcsAffected.size(); j++) {
                int rcID = rcsAffected[j];
                int target = targetNodes[rcID];
                double prevGain, currGain;
                
                // if it is the first time this RC is covered by a seed pair -- below the marginal gain update operations for this case
                if (!isCovered[rcID]) {
                    isCovered[rcID] = true;
                    rcMinSeen[rcID] = std::min(nodeLeanings[target], itemLeanings[bestItemID]);
                    rcMaxSeen[rcID] = std::max(nodeLeanings[target], itemLeanings[bestItemID]);
                    
                    // compute the amount of change in marginal gain for each item contain in this RC set
                    for(int itemID = 0; itemID < nrItems; itemID++) {
                        
                        prevGain = std::abs(nodeLeanings[target] - itemLeanings[itemID]); // mg to empty set
                        if (itemLeanings[itemID] >= rcMinSeen[rcID] &&  itemLeanings[itemID] <= rcMaxSeen[rcID]) {
                            currGain = 0;
                        }
                        else if(itemLeanings[itemID] < rcMinSeen[rcID]) {
                            currGain = rcMinSeen[rcID] - itemLeanings[itemID];
                        }
                        else if(itemLeanings[itemID] > rcMaxSeen[rcID]) {
                            currGain = itemLeanings[itemID] - rcMaxSeen[rcID];
                        }
                        vector<int> rrTemp = rcList->at(itemID)->hyperGT[rcID];
                        for (int z = 0; z < rrTemp.size(); z++) {
                            int pairID = rrTemp[z] + n * itemID;
                            if(pairMark[pairID]) {
                                double tempD = divExpScore[pairID];
                                divExpScore[pairID] = tempD - prevGain + currGain;
                                if(divExpScore[pairID] < 0) {
                                    divExpScore[pairID] = 0;
                                }
                            }//
                        }
                    }
                }
                
                // if it was already covered before
                else {
                    if(itemLeanings[bestItemID] >= rcMinSeen[rcID] && itemLeanings[bestItemID] <= rcMaxSeen[rcID])
                        continue;
                    
                    double prevMin = rcMinSeen[rcID];
                    double prevMax = rcMaxSeen[rcID];
                    rcMinSeen[rcID] = std::min(prevMin, itemLeanings[bestItemID]);
                    rcMaxSeen[rcID] = std::max(prevMax, itemLeanings[bestItemID]);
                    
                    for(int itemID = 0; itemID < nrItems; itemID++) {
                        if (itemLeanings[itemID] >= prevMin &&  itemLeanings[itemID] <= prevMax) {
                            prevGain = 0;
                        }
                        else if(itemLeanings[itemID] < prevMin) {
                            prevGain = prevMin - itemLeanings[itemID];
                        }
                        else if(itemLeanings[itemID] > prevMax) {
                            prevGain = itemLeanings[itemID] - prevMax;
                        }
                        
                        if (itemLeanings[itemID] >= rcMinSeen[rcID] &&  itemLeanings[itemID] <= rcMaxSeen[rcID]) {
                            currGain = 0;
                        }
                        else if(itemLeanings[itemID] < rcMinSeen[rcID]) {
                            currGain = rcMinSeen[rcID] - itemLeanings[itemID];
                        }
                        else if(itemLeanings[itemID] > rcMaxSeen[rcID]) {
                            currGain = itemLeanings[itemID] - rcMaxSeen[rcID];
                        }
                        
                        vector<int> rrTemp = rcList->at(itemID)->hyperGT[rcID];
                        for (int z = 0; z < rrTemp.size(); z++) {
                            int pairID = rrTemp[z] + n * itemID;
                            if(pairMark[pairID]) {
                                double tempD = divExpScore[pairID];
                                divExpScore[pairID] = tempD - prevGain + currGain;
                                if(divExpScore[pairID] < 0) {
                                    divExpScore[pairID] = 0;
                                }
                            }
                        }
                    }
                } // end-if this RC was already covered before
                
            }// end of for loop for affected
        } // end of while k loop
        
        clock_t end = clock();
        totalDuration = (double(end - begin) / CLOCKS_PER_SEC) + duration_common; // in seconds
        totalMemory = getCurrentMemoryUsage(); // in MB
        
        cout << "total degree-far time taken (in seconds) " << totalDuration << " total memory (in mb) " << totalMemory << endl;
        
        string extraFName0 = outFolderName + OS_SEP + "assignment_degreeFar.txt";
        ofstream extraStream0;
        if(extraStream0.is_open())
            extraStream0.close();
        
        extraStream0.open(extraFName0.c_str());
        
        if (extraStream0.is_open() == false) {
            cout << "Can't open file " << extraFName0  << " for writing" << endl;
            exit(1);
        }
        
        extraStream0 << "nodeID" << " " << "itemID" << " " <<  "mgScore" << " " << "cumScore" << " " << "duration(sec)" << " " << "memory(mb)" << endl;
        double totMGs = 0;
        for (int i = 0; i < seedSet.size(); i++) {
            int node = (int) seedSet[i] % n;
            int item = (int) seedSet[i] / n;
            totMGs += seedScores[i];
            extraStream0 << node << " " << item << " " << seedScores[i] << " " << totMGs << " " << totalDuration << " " << totalMemory << endl;
        }
        cout << "total degree-far div. exp. score " << totMGs << endl;
        extraStream0.close();
        
        
        // compute for each seed pair its prob. of reaching target nodes
        
        string extraFName1 = outFolderName + OS_SEP + "pairReachProbs_degreeFar.txt";
        ofstream extraStream1;
        
        if(extraStream1.is_open())
            extraStream1.close();
        
        extraStream1.open(extraFName1.c_str());
        
        if (extraStream1.is_open() == false) {
            cout << "Can't open file " << extraFName1  << " for writing" << endl;
            exit(1);
        }
        
        extraStream1 << "nodeID" << " " << "itemID" << " " <<  "targetNodeID" << " " << "reachProb" << endl;
        
        extraStream1 << std::fixed;
        extraStream1.precision(8);
        
        vector< vector<double> > pairNodeReachProbs;
        
        for (int i = 0; i < k; i++) {
            pairNodeReachProbs.push_back(std::vector<double>(n,0.0));
        }
        
        for (int i = 0; i < k; i++) {
            int seedPairID = seedSet[i];
            vector<int> rcsCovered = hyperG[seedPairID]; // ids of RC sets that contain this best pair
            for (int j = 0; j < rcsCovered.size(); j++) {
                int rcID = rcsCovered[j];
                int target = targetNodes[rcID];
                double tempVal = pairNodeReachProbs[i][target];
                pairNodeReachProbs[i][target] = tempVal + 1.0;
            }
        }
        
        for (int i = 0; i < k; i++) {
            int seedPairID = seedSet[i];
            int nodeID = (int) seedPairID % n;
            int itemID = (int) seedPairID / n;
            for (int j = 0; j < n; j++) {
                double tempVal = pairNodeReachProbs[i][j];
                pairNodeReachProbs[i][j] = (double) n * tempVal / (double) rcSampleSize;
                extraStream1 << nodeID << " " << itemID << " " <<  j << " " << pairNodeReachProbs[i][j] << endl;
            }
        }
        
        extraStream1.close();
        
        // compute fv(seedset) for each target node v
        string extraFName2 = outFolderName + OS_SEP + "nodeExpLevels_degreeFar.txt";
        ofstream extraStream2;
        
        if(extraStream2.is_open())
            extraStream2.close();
        
        extraStream2.open(extraFName2.c_str());
        
        if (!extraStream2.is_open()) {
            cout << "Can't open file " << extraFName2  << " for writing" << endl;
            exit(1);
        }
        
        extraStream2 << "targetNodeID" << " " << "exposure_level" << endl;
        vector<double> nodeExpLevels(n,0.0);
        
        for (int rcID = 0; rcID < rcSampleSize; rcID++) {
            int target = targetNodes[rcID];
            if(isCovered[rcID]) {
                nodeExpLevels[target] += (rcMaxSeen[rcID] - rcMinSeen[rcID]);
            }
        }
        
        //            double kontrolAw = 0;
        for (int i = 0 ; i < n ; i++) {
            double tempVal = nodeExpLevels[i];
            nodeExpLevels[i] = (double) n * tempVal / (double) rcSampleSize;
            //                kontrolAw += nodeExpLevels[i];
            extraStream2 << i << " " << nodeExpLevels[i] << endl;
        }
        
        extraStream2.close();
        //            cout << "kontrol " << kontrolAw << endl;
        
        //        } // end of extra results
        
        
        
        return totalExpScore / (double) rcSampleSize;
    }

    double allocator::degreeWeighted(int64 rcSampleSize, int k) {
        
        clock_t begin = clock();
        priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond>heap;
        std::vector<double> divExpScore(nrPairs, 0);

        std::vector<double> rcMinSeen(rcSampleSize,-5.0);
        std::vector<double> rcMaxSeen(rcSampleSize, 5.0);
        
        vector<bool> isCovered(rcSampleSize, false);
        vector<bool> pairMark(nrPairs, true);
        
        double totalExpScore = 0;
        int bestPairID, bestNodeID, bestItemID;
        
        for (int i = 0; i < n; i++) {
            for (int itemID = 0; itemID < nrItems; itemID++) {
                int prID = i + n * itemID;
                double val = ((double) nodeDegree[i]) * (std::abs(nodeLeanings[i] - itemLeanings[itemID]));
                std::pair<int, double>pairVal(std::make_pair(prID, val));
                heap.push(pairVal);
                divExpScore[prID] = hyper_degree[prID];
            }
        }
        
        
        seedSet.clear();
        seedScores.clear();
        
        
        while ((int)seedSet.size() < k) {
            
            pair<int, double> pairVal = heap.top();
            heap.pop();
            bestPairID = pairVal.first;
            bestNodeID = (int) bestPairID % n;

            
            bestItemID = (int) bestPairID / n;
            totalExpScore += divExpScore[bestPairID];
            seedSet.push_back(bestPairID);
            seedScores.push_back(n * divExpScore[bestPairID] / (double) rcSampleSize);
            pairMark[bestPairID] = false;

            // update the hyper_degree of the pairs whose marginal gain change due to the best pair selection
            vector<int> rcsAffected = hyperG[bestPairID]; // ids of RC sets that contain this best pair
            for (int j = 0; j < rcsAffected.size(); j++) {
                int rcID = rcsAffected[j];
                int target = targetNodes[rcID];
                double prevGain, currGain;
                
                // if it is the first time this RC is covered by a seed pair -- below the marginal gain update operations for this case
                if (!isCovered[rcID]) {
                    isCovered[rcID] = true;
                    rcMinSeen[rcID] = std::min(nodeLeanings[target], itemLeanings[bestItemID]);
                    rcMaxSeen[rcID] = std::max(nodeLeanings[target], itemLeanings[bestItemID]);
                    
                    // compute the amount of change in marginal gain for each item contain in this RC set
                    for(int itemID = 0; itemID < nrItems; itemID++) {
                        
                        prevGain = std::abs(nodeLeanings[target] - itemLeanings[itemID]); // mg to empty set
                        if (itemLeanings[itemID] >= rcMinSeen[rcID] &&  itemLeanings[itemID] <= rcMaxSeen[rcID]) {
                            currGain = 0;
                        }
                        else if(itemLeanings[itemID] < rcMinSeen[rcID]) {
                            currGain = rcMinSeen[rcID] - itemLeanings[itemID];
                        }
                        else if(itemLeanings[itemID] > rcMaxSeen[rcID]) {
                            currGain = itemLeanings[itemID] - rcMaxSeen[rcID];
                        }
                        vector<int> rrTemp = rcList->at(itemID)->hyperGT[rcID];
                        for (int z = 0; z < rrTemp.size(); z++) {
                            int pairID = rrTemp[z] + n * itemID;
                            if(pairMark[pairID]) {
                                double tempD = divExpScore[pairID];
                                divExpScore[pairID] = tempD - prevGain + currGain;
                                if(divExpScore[pairID] < 0) {
                                    divExpScore[pairID] = 0;
                                }
                            }//
                        }
                    }
                }
                
                // if it was already covered before
                else {
                    if(itemLeanings[bestItemID] >= rcMinSeen[rcID] && itemLeanings[bestItemID] <= rcMaxSeen[rcID])
                        continue;
                    
                    double prevMin = rcMinSeen[rcID];
                    double prevMax = rcMaxSeen[rcID];
                    rcMinSeen[rcID] = std::min(prevMin, itemLeanings[bestItemID]);
                    rcMaxSeen[rcID] = std::max(prevMax, itemLeanings[bestItemID]);
                    
                    for(int itemID = 0; itemID < nrItems; itemID++) {
                        if (itemLeanings[itemID] >= prevMin &&  itemLeanings[itemID] <= prevMax) {
                            prevGain = 0;
                        }
                        else if(itemLeanings[itemID] < prevMin) {
                            prevGain = prevMin - itemLeanings[itemID];
                        }
                        else if(itemLeanings[itemID] > prevMax) {
                            prevGain = itemLeanings[itemID] - prevMax;
                        }
                        
                        if (itemLeanings[itemID] >= rcMinSeen[rcID] &&  itemLeanings[itemID] <= rcMaxSeen[rcID]) {
                            currGain = 0;
                        }
                        else if(itemLeanings[itemID] < rcMinSeen[rcID]) {
                            currGain = rcMinSeen[rcID] - itemLeanings[itemID];
                        }
                        else if(itemLeanings[itemID] > rcMaxSeen[rcID]) {
                            currGain = itemLeanings[itemID] - rcMaxSeen[rcID];
                        }
                        
                        vector<int> rrTemp = rcList->at(itemID)->hyperGT[rcID];
                        for (int z = 0; z < rrTemp.size(); z++) {
                            int pairID = rrTemp[z] + n * itemID;
                            if(pairMark[pairID]) {
                                double tempD = divExpScore[pairID];
                                divExpScore[pairID] = tempD - prevGain + currGain;
                                if(divExpScore[pairID] < 0) {
                                    divExpScore[pairID] = 0;
                                }
                            }
                        }
                    }
                } // end-if this RC was already covered before
                
            }// end of for loop for affected
        } // end of while k loop
        
        
        clock_t end = clock();
        totalDuration = (double(end - begin) / CLOCKS_PER_SEC) + duration_common; // in seconds
        totalMemory = getCurrentMemoryUsage(); // in MB
        
        cout << "total degree-weighted time taken (in seconds) " << totalDuration << " total memory (in mb) " << totalMemory << endl;
        
        string extraFName0 = outFolderName + OS_SEP + "assignment_degreeWeighted.txt";
        ofstream extraStream0;
        if(extraStream0.is_open())
            extraStream0.close();
        
        extraStream0.open(extraFName0.c_str());
        
        if (!extraStream0.is_open()) {
            cout << "Can't open file " << extraFName0  << " for writing" << endl;
            exit(1);
        }
        
        extraStream0 << "nodeID" << " " << "itemID" << " " <<  "mgScore" << " " << "cumScore" << " " << "duration(sec)" << " " << "memory(mb)" << endl;
        double totMGs = 0;
        for (int i = 0; i < seedSet.size(); i++) {
            int node = (int) seedSet[i] % n;
            int item = (int) seedSet[i] / n;
            totMGs += seedScores[i];
            extraStream0 << node << " " << item << " " << seedScores[i] << " " << totMGs << " " << totalDuration << " " << totalMemory << endl;
        }
        cout << "total degree-weighted div. exp. score " << totMGs << endl;
        
        extraStream0.close();
        
        
        
        // compute for each seed pair its prob. of reaching target nodes
        
        string extraFName1 = outFolderName + OS_SEP + "pairReachProbs_degreeWeighted.txt";
        ofstream extraStream1;
        
        if(extraStream1.is_open())
            extraStream1.close();
        
        extraStream1.open(extraFName1.c_str());
        
        if (extraStream1.is_open() == false) {
            cout << "Can't open file " << extraFName1  << " for writing" << endl;
            exit(1);
        }
        
        extraStream1 << "nodeID" << " " << "itemID" << " " <<  "targetNodeID" << " " << "reachProb" << endl;
        
        extraStream1 << std::fixed;
        extraStream1.precision(8);
        
        vector< vector<double> > pairNodeReachProbs;
        
        for (int i = 0; i < k; i++) {
            pairNodeReachProbs.push_back(std::vector<double>(n,0.0));
        }
        
        for (int i = 0; i < k; i++) {
            int seedPairID = seedSet[i];
            vector<int> rcsCovered = hyperG[seedPairID]; // ids of RC sets that contain this best pair
            for (int j = 0; j < rcsCovered.size(); j++) {
                int rcID = rcsCovered[j];
                int target = targetNodes[rcID];
                double tempVal = pairNodeReachProbs[i][target];
                pairNodeReachProbs[i][target] = tempVal + 1.0;
            }
        }
        
        for (int i = 0; i < k; i++) {
            int seedPairID = seedSet[i];
            int nodeID = (int) seedPairID % n;
            int itemID = (int) seedPairID / n;
            for (int j = 0; j < n; j++) {
                double tempVal = pairNodeReachProbs[i][j];
                pairNodeReachProbs[i][j] = (double) n * tempVal / (double) rcSampleSize;
                extraStream1 << nodeID << " " << itemID << " " <<  j << " " << pairNodeReachProbs[i][j] << endl;
            }
        }
        
        extraStream1.close();
        
        // compute fv(seedset) for each target node v
        string extraFName2 = outFolderName + OS_SEP + "nodeExpLevels_degreeWeighted.txt";
        ofstream extraStream2;
        
        if(extraStream2.is_open())
            extraStream2.close();
        
        extraStream2.open(extraFName2.c_str());

        if (extraStream2.is_open()) {

            extraStream2 << "targetNodeID" << " " << "exposure_level" << endl;
            vector<double> nodeExpLevels(n, 0.0);

            for (int rcID = 0; rcID < rcSampleSize; rcID++) {
                int target = targetNodes[rcID];
                if (isCovered[rcID]) {
                    nodeExpLevels[target] += (rcMaxSeen[rcID] - rcMinSeen[rcID]);
                }
            }

            //            double kontrolAw = 0;
            for (int i = 0; i < n; i++) {
                double tempVal = nodeExpLevels[i];
                nodeExpLevels[i] = (double) n * tempVal / (double) rcSampleSize;
                //                kontrolAw += nodeExpLevels[i];
                extraStream2 << i << " " << nodeExpLevels[i] << endl;
            }

            extraStream2.close();
            //            cout << "kontrol " << kontrolAw << endl;

            //        } // end of extra results



            return totalExpScore / (double) rcSampleSize;
        } else {
            cout << "Can't open file " << extraFName2 << " for writing" << endl;
            exit(1);
        }
    }

    void allocator::readTICGraph() {
        
        string probGraphFile = opt->getValue("probGraphFile");
        cout << "Reading file " << probGraphFile << endl;
        ifstream myfile (probGraphFile.c_str(), ios::in);
        
        double *itemProbs; //item-specific influence probabilities
        
        int nrEdges = 0;
        set<int> nodes; // for control
        
        if (myfile.is_open()) {
            while (! myfile.eof() )	{
                std::string line;
                getline (myfile,line);
                if (line.empty()) continue;
                
                std::string::size_type pos = line.find_first_of(delim);
                int prevpos = 0;
                
                //first user
                string str = line.substr(prevpos, pos-prevpos);
                int u1 = strToInt(str);
                
                //second user
                prevpos = line.find_first_not_of(delim, pos);
                pos = line.find_first_of(delim, prevpos);
                int u2 = strToInt(line.substr(prevpos, pos-prevpos));
                
                if (u1 == u2)
                    continue;
                
                nodeDegree[u1] = nodeDegree[u1] + 1.0;
                
                nrEdges++;
                
                graphT[u2].push_back(u1); //insert to the transposed graph
                
                // for control
                nodes.insert(u1);
                nodes.insert(u2);
                
                prevpos = line.find_first_not_of(delim, pos);
                
                str = line.substr(prevpos);
                itemProbs = new double[nrItems];
                stringTokenizer(str, itemProbs, nrItems, delim);
                
                
                for(int i = 0; i < nrItems; i++) {
                    rcList->at(i)->probT[u2].push_back(itemProbs[i]);
                }
                
            }
            
            // verify input
            if (n != (int) nodes.size()) {
                cout << "problem: nr of nodes read is not equal to input n, exiting..." << endl;
                exit(1);
            }
            
            if (m != (int) nrEdges) {
                cout << "problem: nr of edges read is not equal to input m, exiting..." << endl;
                exit(1);
            }
            
//            cout << "node ids kontrol " << endl;
//            for (set<int>::iterator it = nodes.begin(); it != nodes.end(); it++) {
//                cout << (*it) << " ";
//            }
//            cout << endl;
            
            myfile.close();
        }
        
        else
            cout << "Can't open input graph file " << probGraphFile << endl;
        
        cout << "graph import complete " << endl;
        cout << "total number of nodes " << nodes.size() << endl;
        cout << "total number of edges " << nrEdges << endl;
    }
    
    void allocator::readItemLeaningsFile() {
        // temporary, to delete later
        itemLeanings.push_back(strToDouble("0"));
        itemLeanings.push_back(strToDouble("1"));
    }
    
    void allocator::readNodeLeaningsFile() {
        // temporary, to delete later
        int nodeIndex = 0;
        if(++nodeIndex != n){
            nodeLeanings.push_back(strToDouble("0.0"));
        }
    }
    
    allocator::~allocator() {
        cout << "assignments complete! " << endl;
    }
    
    void allocator::arrangeOutputFiles() {
        
                string command = string("mkdir -p ") + outFolderName ;
        
                system(command.c_str());
        
        
        string masterFileName = "assignment_tdem.txt";
        outMasterName = outFolderName + OS_SEP + masterFileName;
//        outMasterName = masterFileName;
        
        if(outMasterStream.is_open())
            outMasterStream.close();
        
        outMasterStream.open(outMasterName.c_str());
        
        outMasterStream << "nodeID" << " " << "itemID" << " " <<  "mgScore" << " " << "cumScore" << " " << "duration(sec)" << " " << "memory(mb)" << endl;
        
        if (outMasterStream.is_open() == false) {
            cout << "Can't open file " << outMasterName  << " for writing" << endl;
            exit(1);
        }
        
        // memory icin
        command = string("mkdir -p temp") ;
        system(command.c_str());
    }
    
    //    void allocator::writeInMasterOutputFile(int nodeID, int itemID) {
    //        // seed-node item mgScore totScore runTime(sec) memory(mb)
    //        outMasterStream << nodeID << " " << itemID << endl;
    //    }
    
    
    void allocator::writeInMasterOutputFile(int nodeRed, int nodeBlue, double mgScore, double totScore, float duration, float memory) {
        // seed-node item mgScore totScore runTime(sec) memory(mb)
        outMasterStream << nodeRed << " " << nodeBlue << " " <<  mgScore << " " << totScore << " " << duration << " " << memory << endl;
    }

    
}
