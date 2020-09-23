#include "allocator.h"
#include "anyoption.h"
#include "iheap.h"
#include <ctime>
#include <numeric>

namespace _Cide{
    
    allocator::allocator(AnyOption* opt1, string inputname, string comparefilename) {
        
        opt = opt1;
        probGraphFile = inputname;
        //delim = " \t";
        delim = " ";

        compareFile = comparefilename;
        compareFileShort = compareFile;
        for (int j = 0; j < 4; ++j) {
           compareFileShort.pop_back();
        }

        readGraphNodes(); // assign n and m
        readCompareNodes(); // get red and blue nodes as k_r and k_b

        sfmt_init_gen_rand(&sfmtSeed , 95082);
        
        //n = strToInt(opt->getValue("n"));
        //m = strToInt(opt->getValue("m"));

        //k = strToInt(opt->getValue("k"));
        //k_r = strToInt(opt->getValue("kr"));
        //k_b = strToInt(opt->getValue("kb"));
        tao = ceil((double)k_b/(double)k_r);
        P = Permutation(n, k_r*(tao + 1)) / (double) (pow(fact(tao), k_r) * fact(k_r));


        nrItems = strToInt(opt->getValue("nrItems"));
        epsilon = strToDouble(opt->getValue("epsilon"));
        ell = strToDouble(opt->getValue("ell"));
        outFolderName = opt->getValue("outputFolder");
    
        nrPairs = n * nrItems; // multiple items

        nodeDegree.resize(n);

        
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
            //for (int j = 0; j < n; j++)
                //ig->probT.push_back(std::vector< double>());
                //ig->probT.emplace_back(std::vector< double>());
            rcList->push_back(ig);
        }

        //nodeDegree = std::vector<double>(n,0.0);
        
        cout << "nr items " << nrItems << endl;
        //cout << "assignment size " << k << endl;

        readTICGraph();

        cout << "successfully read everything " << endl;

        arrangeOutputFiles();
        cout << "starting assignment for tdem" << endl;
        tdem();
        
// comment out below to run the baselines
        degreeVersionOne(prevSize);
        degreeVersionTwo(prevSize);
        compareGivenNodes(prevSize);
    }


    double allocator::tdem() {
        
        cout << "computation for lower bounding OPT started " << endl;

        clock_t common_begin = clock();
        double lb = lowerBoundOPT();
        cout << "lower bound identified: " << lb << endl;
        
        theta = (2 + 2/3.0 * epsilon) * n * (log(n) + log(2) +  log(P)) / (epsilon * epsilon * lb);
        
        cout << "final sample size " << theta << endl;
        generateRCSets(theta);
        
        clock_t common_end = clock();
        duration_common = double(common_end - common_begin) / CLOCKS_PER_SEC;

        double greedySolution = n * rcGreedy(prevSize, true);

        cout <<"Greedy Solution is: "<< greedySolution << endl;

        return greedySolution; 
    }

    double allocator::lowerBoundOPT() {

        double epsilon_1 = epsilon;

        for (int x = 1; x < log2(n); x++) {

            int64 theta_x = (2+2.0/3.0 * epsilon_1)* (ell * log(n) + log(P) + log(log2(n))) * pow(2.0,x) / (epsilon_1 * epsilon_1);
            generateRCSets(theta_x);
            double ept = rcGreedy(theta_x, false);

            cout << "ept is: " << ept << endl;

            if (ept > ((1+epsilon_1) * 2.0 / pow(2.0, x))) {
                double lowerBound = ept * n / (1 + epsilon_1);
                return lowerBound;
            }
        }
        cout << "returning naive lower bound  " << endl;
        double naive = 1.0;
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
        cout << "current size is " << newSize <<  " size of hyPerGTPairs is " << tmpsize << endl;

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
                hyperG[y].push_back(rcID); // careful about the repetition
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
        //cout << "newSize is " << newSize << endl;
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


    double allocator::rcGreedy(int64 rcSampleSize, bool extraResults) {
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

        double totalExpScore = 0;
        int bestPairID, bestRedNodeID, bestBlueNodeID;
        std::vector<int> isValidRed = std::vector<int>(n,tao);
        std::vector<int> isValidBlue = std::vector<int>(n,1);
        std::vector<bool> isCovered = std::vector<bool>(rcSampleSize, false);


        seedSet.clear();
        seedScores.clear(); // seedScores should be cleared as well.
        std::set<int> rednodes;
        int flag_red = 0;
        while ((int)seedSet.size() < k_b) {
            // here I use the vector to store the results, but can use the priority queue.
            //bestPairID = std::max_element(ExpScore.begin(), ExpScore.end()) - ExpScore.begin();
            bestPairID = std::distance(ExpScore.begin(), std::max_element(ExpScore.begin(), ExpScore.end()));
            //cout << "bestPairID is " << bestPairID   <<endl;
            bestRedNodeID = (int) bestPairID / n;
            rednodes.insert(bestRedNodeID);
            bestBlueNodeID = (int) bestPairID % n;

            cout << "sum of Expscore is " << accumulate(ExpScore.begin(), ExpScore.end(), 0) << endl;


            // if can not select, then just pop out
            if (isValidRed[bestRedNodeID] == 0 || isValidBlue[bestBlueNodeID] == 0) {
                ExpScore[bestPairID] = 0;
                continue;
            }

            // push
            seedSet.push_back(bestPairID); // store
            //cout << "ExpScore is " << ExpScore[bestPairID] << "rc sample size is " <<rcSampleSize << "n is " << n << endl;
            seedScores.push_back(ExpScore[bestPairID]); // store score, dividing later

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

            // if red nodes reaches k_r
            if(rednodes.size() == k_r && flag_red == 0){
                for (int i = 0; i < n; ++i) {
                    if(rednodes.count(i) == 0){
                        isValidRed[i] = 0;
                        for (int j = 0; j < n; ++j) {
                           ExpScore[i*n + j] = 0;
                        }
                    }
                }
                flag_red = 1;
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

        int covercount = 0;
        for (int l = 0; l < isCovered.size(); ++l) {
            if(isCovered[l])
                covercount ++;
        }

        cout << "totalExpScore is: " << totalExpScore << " cover is " << covercount << " sample size is " << rcSampleSize << endl;


        ///////////// produce the additional results here
        ////////////  Remember should write the f function results not the g function result.
        if(extraResults) {

            clock_t end = clock();
            totalDuration = double(end - begin) / CLOCKS_PER_SEC + duration_common; // in seconds
            totalMemory = getCurrentMemoryUsage(); // in MB
            
            //cout << "total tdem time taken (in seconds) " << totalDuration << " total memory (in mb) " << totalMemory << endl;
            
            // compute for each seed pair its prob. of reaching target nodes
            
            // compute f(seedset) for each target node v, should consider the other pairs.
            string extraFName0 = outFolderName + OS_SEP + compareFileShort + OS_SEP + "Greedy.txt";
            ofstream extraStream0;

            if(extraStream0.is_open())
                extraStream0.close();

            extraStream0.open(extraFName0.c_str());

            if (! extraStream0.is_open()) {
                cout << "Can't open file " << extraFName0  << " for writing" << endl;
                exit(1);
            }

            vector<double> nodeSampleExpLevels(rcSampleSize,0.0);
            vector<double> nodeExpLevels(n,0.0);
            std::vector<int> rrRed;
            std::vector<int> rrBlue;

            std::vector<int> redseedSet;
            std::vector<int> blueseedSet;
            redseedSet.clear();
            blueseedSet.clear();

            float tmptotalscore = 0;
            for (int i = 0; i < k_b; ++i) {
                int SeedPairID = seedSet[i];
                tmptotalscore += seedScores[i];
                //writeInMasterOutputFile(SeedPairID/n, SeedPairID%n, seedScores[i], tmptotalscore, totalDuration, totalMemory);
                for (int id : hyperG[SeedPairID / n]) {
                    rrRed.push_back(id);
                }
                redseedSet.push_back(SeedPairID/n);
                for (int id : hyperG[SeedPairID % n + n]) {
                    rrBlue.push_back(id);
                }
                blueseedSet.push_back(SeedPairID%n);
            }

            extraStream0 << "red nodes: ";
            for (int k = 0; k < redseedSet.size(); ++k) {
                extraStream0 << redseedSet[k] << " ";
            }
            extraStream0 << "\n";

            extraStream0 << "blue nodes: ";
            for (int k = 0; k < blueseedSet.size(); ++k) {
                extraStream0 << blueseedSet[k] << " ";
            }
            extraStream0 << "\n";

            std::sort(rrRed.begin(), rrRed.end());
            std::sort(rrBlue.begin(), rrBlue.end());

            std::vector<int>::iterator it;
            std::vector<int> rrRedBlue(rrRed.size() + rrBlue.size(), 0);

            it = std::set_intersection(rrRed.begin(), rrRed.end(), rrBlue.begin(), rrBlue.end(), rrRedBlue.begin());
            rrRedBlue.resize(it - rrRedBlue.begin());

            float totalscore = (float) rrRedBlue.size() / (float) rcSampleSize;

            extraStream0 << "Time" << " " << "Memory (in mb)"<< endl;
            extraStream0 << totalDuration << " " << totalMemory <<endl;

            extraStream0 <<  "Score is " << totalscore*n << endl;
            extraStream0.close();
//
            cout << "Total score of Greedy is " << totalscore*n << endl;

            writeInMasterOutputFile("Greedy", totalDuration,  totalMemory, totalscore*n);

        } // end of extra results
        
        return totalExpScore / (double)rcSampleSize;
    }


    double allocator::degreeVersionOne(int64 rcSampleSize) {

        string extraFName0 = outFolderName + OS_SEP + compareFileShort + OS_SEP + "degreeVersionOne.txt";
        ofstream extraStream0;
        if(extraStream0.is_open())
            extraStream0.close();

        extraStream0.open(extraFName0.c_str());

        if (!extraStream0.is_open()) {
            cout << "Can't open file " << extraFName0  << " for writing" << endl;
            exit(1);
        }


        clock_t begin = clock();

        priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond> heap;

        float totalscore = 0;
        int bestNodeID;
//
//
        for (int i = 0; i < n; i++) {
            double val = (double) nodeDegree[i];
            std::pair<int, double> pairVal(std::make_pair(i, val));
            heap.push(pairVal);
        }

        seedSet.clear();
        seedScores.clear();
        vector<double> nodeSampleExpLevels(rcSampleSize, 0.0);
        vector<double> nodeExpLevels(n, 0.0);
        std::vector<int> rrRed;
        std::vector<int> rrBlue;

        extraStream0 << "red nodes: ";

        while ((int) seedSet.size() < k_r) {
            pair<int, double> pairVal = heap.top();
            heap.pop();
            bestNodeID = pairVal.first;
            for (int id : hyperG[bestNodeID]) {
                rrRed.push_back(id);
            }
            seedSet.push_back(bestNodeID);
            extraStream0 << bestNodeID << " ";
        }
        extraStream0 << "\n";

        extraStream0 << "blue nodes: ";

        while ((int) seedSet.size() < k_r + k_b) {
            pair<int, double> pairVal = heap.top();
            heap.pop();
            bestNodeID = pairVal.first;
            for (int l = 0; l < hyperG[bestNodeID].size(); ++l) {
                int id = hyperG[bestNodeID][l];
                rrBlue.push_back(id);
            }
            seedSet.push_back(bestNodeID);
            extraStream0 << bestNodeID << " ";
        }
        extraStream0 << "\n";

        std::sort(rrRed.begin(), rrRed.end());
        std::sort(rrBlue.begin(), rrBlue.end());

        std::vector<int>::iterator it;
        std::vector<int> rrRedBlue(rrRed.size() + rrBlue.size(), 0);

        it = std::set_intersection(rrRed.begin(), rrRed.end(), rrBlue.begin(), rrBlue.end(), rrRedBlue.begin());
        rrRedBlue.resize(it - rrRedBlue.begin());

        totalscore = (float) rrRedBlue.size()/ (float) rcSampleSize;

        // end of extra results
        clock_t end = clock();
        totalDuration = double(end - begin) / CLOCKS_PER_SEC + duration_common; // in seconds
        totalMemory = getCurrentMemoryUsage(); // in MB

        extraStream0 << "Time" << " " << "Memory (in mb)"<< endl;
        extraStream0 << totalDuration << " " << totalMemory <<endl;

        cout << "Total score of Degree one is" << " " << totalscore*n << endl;
        extraStream0 << "score is" << " " << totalscore*n << endl;

        extraStream0.close();

        writeInMasterOutputFile("DegreeVersionOne", totalDuration,  totalMemory, totalscore*n);

        return totalscore;
    }


    double allocator::degreeVersionTwo(int64 rcSampleSize) {

        string extraFName0 = outFolderName + OS_SEP + compareFileShort + OS_SEP + "degreeVersionTwo.txt";

        ofstream extraStream0;
        if(extraStream0.is_open())
            extraStream0.close();

        extraStream0.open(extraFName0.c_str());

        if (!extraStream0.is_open()) {
            cout << "Can't open file " << extraFName0  << " for writing" << endl;
            exit(1);
        }

        clock_t begin = clock();
        priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond> heap;

        float totalscore = 0;
        int bestNodeID;
//
        for (int i = 0; i < n; i++) {
            double val = (double) nodeDegree[i];
            std::pair<int, double> pairVal(std::make_pair(i, val));
            heap.push(pairVal);
        }

        std::vector<int> redseedSet, blueseedSet;
        redseedSet.clear();
        blueseedSet.clear();
        seedSet.clear();
        seedScores.clear();
        vector<double> nodeSampleExpLevels(rcSampleSize, 0.0);
        vector<double> nodeExpLevels(n, 0.0);
        std::vector<int> rrRed;
        std::vector<int> rrBlue;

        //rrRed.reserve(n*k_r);
        //rrBlue.reserve(n*k_b);

        while ((int) seedSet.size() < k_r + k_r) {
            pair<int, double> pairVal = heap.top();
            heap.pop();
            bestNodeID = pairVal.first;
            for (int j = 0; j < hyperG[bestNodeID].size(); ++j) {
                int id = hyperG[bestNodeID][j];
                rrRed.push_back(id);
            }

            seedSet.push_back(bestNodeID);
            redseedSet.push_back(bestNodeID);

            pairVal = heap.top();
            heap.pop();
            bestNodeID = pairVal.first;
            for (int j = 0; j < hyperG[bestNodeID].size(); ++j) {
                int id = hyperG[bestNodeID][j];
                rrBlue.push_back(id);
            }

            seedSet.push_back(bestNodeID);
            blueseedSet.push_back(bestNodeID);
        }

        while ((int) seedSet.size() < k_r + k_b) {
            pair<int, double> pairVal = heap.top();
            heap.pop();
            bestNodeID = pairVal.first;
            for (int l = 0; l < hyperG[bestNodeID].size(); ++l) {
                int id = hyperG[bestNodeID][l];
                rrBlue.push_back(id);
            }
            seedSet.push_back(bestNodeID);
            blueseedSet.push_back(bestNodeID);
        }

        extraStream0 << "red nodes: ";
        for (int k = 0; k < redseedSet.size(); ++k) {
            extraStream0 << redseedSet[k] << " ";
        }
        extraStream0 << "\n";

        extraStream0 << "blue nodes: ";
        for (int k = 0; k < blueseedSet.size(); ++k) {
            extraStream0 << blueseedSet[k] << " ";
        }
        extraStream0 << "\n";

        std::sort(rrRed.begin(), rrRed.end());
        std::sort(rrBlue.begin(), rrBlue.end());

        std::vector<int>::iterator it;
        std::vector<int> rrRedBlue(rrRed.size() + rrBlue.size(), 0);

        it = std::set_intersection(rrRed.begin(), rrRed.end(), rrBlue.begin(), rrBlue.end(), rrRedBlue.begin());
        rrRedBlue.resize(it - rrRedBlue.begin());

        totalscore = (float) rrRedBlue.size() / (float) rcSampleSize;

        // end of extra results
        clock_t end = clock();
        totalDuration = double(end - begin) / CLOCKS_PER_SEC + duration_common; // in seconds
        totalMemory = getCurrentMemoryUsage(); // in MB

        extraStream0 << "Time" << " " << "Memory (in mb)"<< endl;
        extraStream0 << totalDuration << " " << totalMemory <<endl;

        extraStream0 <<  "Score is " << totalscore*n << endl;

        extraStream0.close();

        cout << "Total score of degree two is " << totalscore*n << endl;

        writeInMasterOutputFile("DegreeVersionTwo", totalDuration,  totalMemory, totalscore*n);
        return totalscore;
    }

    double allocator::compareGivenNodes(int64 rcSampleSize) {

        string extraFName0 = outFolderName + OS_SEP + compareFileShort + OS_SEP + "compareGiven.txt";

        ofstream extraStream0;
        if(extraStream0.is_open())
            extraStream0.close();

        extraStream0.open(extraFName0.c_str());

        if (!extraStream0.is_open()) {
            cout << "Can't open file " << extraFName0  << " for writing" << endl;
            exit(1);
        }

        clock_t begin = clock();

        std::vector<int> rrRed;
        std::vector<int> rrBlue;

        for (int i : comparered) {
            for (int id : hyperG[i]) {
               rrRed.push_back(id);
            }
        }

        for (int i : compareblue) {
            for (int id : hyperG[i]) {
                rrBlue.push_back(id);
            }
        }


        extraStream0 << "red nodes: ";
        for (int k : comparered) {
            extraStream0 << k << " ";
        }
        extraStream0 << "\n";

        extraStream0 << "blue nodes: ";
        for (int k : compareblue) {
            extraStream0 << k << " ";
        }
        extraStream0 << "\n";

        std::sort(rrRed.begin(), rrRed.end());
        std::sort(rrBlue.begin(), rrBlue.end());

        std::vector<int>::iterator it;
        std::vector<int> rrRedBlue(rrRed.size() + rrBlue.size(), 0);

        it = std::set_intersection(rrRed.begin(), rrRed.end(), rrBlue.begin(), rrBlue.end(), rrRedBlue.begin());
        rrRedBlue.resize(it - rrRedBlue.begin());

        float totalscore = (float) rrRedBlue.size() / (float) rcSampleSize;

        // end of extra results
        clock_t end = clock();
        totalDuration = double(end - begin) / CLOCKS_PER_SEC + duration_common; // in seconds
        totalMemory = getCurrentMemoryUsage(); // in MB

        extraStream0 << "Time" << " " << "Memory (in mb)"<< endl;
        extraStream0 << totalDuration << " " << totalMemory <<endl;

        extraStream0 <<  "Score is " << totalscore*n << endl;
        extraStream0.close();

        cout << "Total score of compare is " << totalscore*n << endl;


        writeInMasterOutputFile("Compare", totalDuration,  totalMemory, totalscore*n);

        return totalscore;
    }


    void allocator::readCompareNodes() {

        //string probGraphFile = opt->getValue("probGraphFile");
        ifstream myfile (compareFile.c_str(), ios::in);

        //double *itemProbs; //item-specific influence probabilities

        if (myfile.is_open()) {
            string line1;
            string line2;
            getline(myfile, line1);
            getline(myfile, line2);

            vector<string> nodeids1;
            vector<string> nodeids2;

            split1(line1, nodeids1);
            split1(line2, nodeids2);

            if(nodeids1.size() > nodeids2.size()){
                for (int i = 2; i < nodeids1.size(); ++i) {
                    unsigned int tt = strToInt(nodeids1[i]);
                    compareblue.push_back(tt);
                }
                for (int i = 2; i < nodeids2.size(); ++i) {
                    unsigned int tt = strToInt(nodeids2[i]);
                    comparered.push_back(tt);
                }
            }
            else{
                for (int i = 2; i < nodeids1.size(); ++i) {
                    unsigned int tt = strToInt(nodeids1[i]);
                    comparered.push_back(tt);
                }
                for (int i = 2; i < nodeids2.size(); ++i) {
                    unsigned int tt = strToInt(nodeids2[i]);
                    compareblue.push_back(tt);
                }
            }

            myfile.close();
        }

        k_r = comparered.size();
        k_b = compareblue.size();
        cout << "total number of red nodes " << k_r << endl;
        cout << "total number of blue nodes " << k_b << endl;
    }



    void allocator::readGraphNodes() {

        //string probGraphFile = opt->getValue("probGraphFile");
        ifstream myfile (probGraphFile.c_str(), ios::in);

        //double *itemProbs; //item-specific influence probabilities

        int nrEdges = 0;
        set<int> nodes; // for control

        if (myfile.is_open()) {
            while (! myfile.eof() )	{
                std::string line;
                getline (myfile,line);
                if (line.empty()) continue;

                std::vector<std::string> nodeids;
                split1(line, nodeids);
                unsigned int u1 = strToInt(nodeids[0]);
                unsigned int u2 = strToInt(nodeids[1]);

                if (u1 == u2)
                    continue;

                nrEdges++;
                // for control
                nodes.insert(u1);
                nodes.insert(u2);

            }

            myfile.close();
        }

        n = nodes.size();
        m = nrEdges;
        cout << "total number of nodes " << nodes.size() << endl;
        cout << "total number of edges " << nrEdges << endl;
    }



    void allocator::readTICGraph() {
        
        //string probGraphFile = opt->getValue("probGraphFile");
        cout << "Reading file " << probGraphFile << endl;
        ifstream myfile (probGraphFile.c_str(), ios::in);


        //double *itemProbs; //item-specific influence probabilities

        int nrEdges = 0;
        set<int> nodes; // for control
        
        if (myfile.is_open()) {
            while (! myfile.eof() )	{
                std::string line;
                getline (myfile,line);
                if (line.empty()) continue;

                std::vector<std::string> nodeids;
                split1(line, nodeids);
                unsigned int u1 = strToInt(nodeids[0]);
                unsigned int u2 = strToInt(nodeids[1]);
                
                if (u1 == u2)
                    continue;

                nodeDegree[u1] = nodeDegree[u1] + 1.0;

                nrEdges++;
                
                graphT[u2].push_back(u1); //insert to the transposed graph
                
                // for control
                nodes.insert(u1);
                nodes.insert(u2);

            }


            myfile.close();
        }
        
        else
            cout << "Can't open input graph file " << probGraphFile << endl;
        
        cout << "graph import complete " << endl;
    }

    
    allocator::~allocator() {
        cout << "assignments complete! " << endl;
    }
    
    void allocator::arrangeOutputFiles() {
        
       string command = string("mkdir -p ") + outFolderName ;
        
       system(command.c_str());

        string commands = string("mkdir -p ") + outFolderName + OS_SEP + compareFileShort ; // results in detail
        system(commands.c_str());
        
        string masterFileName = compareFile + "_combinedresults.txt";
        outMasterName = outFolderName + OS_SEP + masterFileName;
//        outMasterName = masterFileName;
        
        if(outMasterStream.is_open())
            outMasterStream.close();
        
        outMasterStream.open(outMasterName.c_str());
        
        outMasterStream << "n" << " " << "m" << " " << "redNodes" << " " << "blueNodes" << " "<< endl;
        outMasterStream << n << " " << m << " " << k_r << " " << k_b << " " << endl;
        outMasterStream << "Algorithm" << " " << "Time" << " " << "Memory(MB)" << " " << "score" << endl;
        
        if (!outMasterStream.is_open()) {
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
    
    
    void allocator::writeInMasterOutputFile(const string& algorithm, float duration, float memory, double score) {
        // seed-node item mgScore totScore runTime(sec) memory(mb)
        outMasterStream << algorithm << " " << duration << " " <<  memory << " " << score <<  endl;
    }


}
