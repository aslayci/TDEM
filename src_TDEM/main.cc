#include <cstdlib>
#include <iostream>
#include "anyoption.h"
#include "allocator.h"


AnyOption* readOptions(int argc, char* argv[]);

using namespace _Cide;

string inputname;
string comparename;

int main(int argc, char* argv[]) {
	
    AnyOption *opt = readOptions(argc,argv);
    _Cide::allocator *alloc = new _Cide::allocator(opt, inputname, comparename);
    delete alloc;
    
}

AnyOption* readOptions(int argc, char* argv[]) {
    
    AnyOption *opt = new AnyOption();
    
    // ignore POSIX style options
    opt->noPOSIX();
    opt->setVerbose(); /* print warnings about unknown options */
    opt->autoUsagePrint(true); /* print usage for bad options */
    
    opt->addUsage("");
    opt->addUsage("Usage: ");
    opt->addUsage("");
    opt->addUsage("-help Prints this help ");
    opt->addUsage(" -c <config_file> Specify config file ");
    opt->addUsage(" -x <input_file> Specify input file name");
    opt->addUsage("");
    
    opt->setOption("probGraphFile");
    opt->setOption("itemLeaningsFile");
    opt->setOption("outputFolder");
    
    //opt->setOption("n");
    //opt->setOption("m");
    opt->setOption("nrItems");
    //opt->setOption("k");
    //opt->setOption("kr");
    //opt->setOption("kb");
//    opt->setOption("attentionConstraint");

    opt->setOption("epsilon");
    opt->setOption("ell");

    opt->setCommandFlag("help");
    opt->setCommandOption("c");
    opt->setCommandOption("x"); // to specify input file
    opt->setCommandOption("y"); // to specify comparing file
    opt->processCommandArgs(argc, argv);


    inputname = opt->getValue("x");
    if (inputname.empty()) {
        cout << "input file not memtioned" << endl;
        opt->printUsage();
        delete opt;
        exit(0);
    }
    cout<<"inputfile name is " << inputname << endl;

    comparename = opt->getValue("y");
    if (comparename.empty()) {
        cout << "compare file not memtioned" << endl;
        opt->printUsage();
        delete opt;
        exit(0);
    }
    cout<<"comparfile name is " << comparename << endl;

    if(opt->getFlag( "help" )) {
        opt->printUsage();
        delete opt;
        exit(0);
    }
    
    const char* configFile = opt->getValue("c");
    if (configFile == NULL) {
        cout << "Config file not mentioned" << endl;
        opt->printUsage();
        delete opt;
        exit(0);
    }
    
    cout << "Config file : " << configFile << endl;
    opt->processFile(configFile);
    opt->processCommandArgs( argc, argv );
    cout << endl << endl;
    return opt;
}
