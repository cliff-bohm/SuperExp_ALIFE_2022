//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#pragma once    // directive to insure that this .h file is only included one time

#include <World/AbstractWorld.h> // AbstractWorld defines all the basic function templates for worlds
#include <string>
#include <memory> // shared_ptr
#include <map>

//using std::shared_ptr;
//using std::string;
//using std::map;
//using std::unordered_map;
//using std::unordered_set;
//using std::to_string;

class SawGenesWorld : public AbstractWorld {

public:
    // parameters for group and brain namespaces
    static std::shared_ptr<ParameterLink<std::string>> modePL;
    static std::shared_ptr<ParameterLink<int>> numGenesPL;
    static std::shared_ptr<ParameterLink<int>> geneSizePL;
    static std::shared_ptr<ParameterLink<std::string>> geneToScorePL;
    static std::shared_ptr<ParameterLink<std::string>> popSizePL;

    std::string mode;
    int numGenes, geneSize;
    std::vector<double> geneToScore;

    std::vector<int> popSizes = { -1 }; // default to use initPop
    std::vector<int> popTimes = { -1 }; // default to not change
    int popSizeIndex = 0;


    // a local variable used for faster access to the ParameterLink value
    
    std::string groupName = "root::";
    std::string genomeName = "root::";
    
    SawGenesWorld(std::shared_ptr<ParametersTable> PT);
    virtual ~SawGenesWorld() = default;

    virtual auto evaluate(std::map<std::string, std::shared_ptr<Group>>& /*groups*/, int /*analyze*/, int /*visualize*/, int /*debug*/) -> void override;

    virtual auto requiredGroups() ->std::unordered_map<std::string, std::unordered_set<std::string>> override;
};

