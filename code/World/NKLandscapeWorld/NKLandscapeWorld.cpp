//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

// Evaluates agents on how many '1's they can output. This is a purely fixed
// task
// that requires to reactivity to stimuli.
// Each correct '1' confers 1.0 point to score, or the decimal output determined
// by 'mode'.

#include "NKLandscapeWorld.h"
#include "../../Utilities/Random.h"
#include <cmath>
#include <iostream>
#include <string>

std::shared_ptr<ParameterLink<int>> NKLandscapeWorld::modePL = Parameters::register_parameter("WORLD_NKLANDSCAPE-mode", 0, "0 = bit outputs before adding, 1 = add outputs");
std::shared_ptr<ParameterLink<int>> NKLandscapeWorld::N_PL = Parameters::register_parameter("WORLD_NKLANDSCAPE-n", 1,"A number that defines the genome size.");
std::shared_ptr<ParameterLink<int>> NKLandscapeWorld::K_PL = Parameters::register_parameter("WORLD_NKLANDSCAPE-k", 1,"A number between 1 and N(genome size)");
std::shared_ptr<ParameterLink<std::string>> NKLandscapeWorld::groupNamePL = Parameters::register_parameter("WORLD_NKLANDSCAPE_NAMES-groupNameSpace",(std::string) "root::","namespace of group to be evaluated");
std::shared_ptr<ParameterLink<std::string>> NKLandscapeWorld::genomeNamePL = Parameters::register_parameter("WORLD_NKLANDSCAPE_NAMES-genomeNameSpace", (std::string) "root::","namespace for parameters used to define brain");
std::shared_ptr<ParameterLink<int>> NKLandscapeWorld::numTraitsPL = Parameters::register_parameter("WORLD_NKLANDSCAPE-numTraits",1, "Number of redundant copies to include in one organizm. They are all independant but all contribute to fitness.");

NKLandscapeWorld::NKLandscapeWorld(std::shared_ptr<ParametersTable> PT_)
    : AbstractWorld(PT_) {
    
    std::cout << "Generating landscape (const seed 1337)" << std::endl;
    std::mt19937 rng(1337);
    landscape = std::vector<std::vector<double>>( N_PL->get(PT), std::vector<double>() );
    for (int i=0; i<N_PL->get(PT); i++){
        for (int j=0; j < std::pow(2, K_PL->get(PT)); j++){
            auto value = Random::getDouble(1.0, rng);
	    landscape[i].push_back(value);
            std::cout << value;
        }
        std::cout << std::endl;
    }
    

    popFileColumns.clear();

    genomeNames = std::unordered_set<std::string>();
    for (int i = 0; i < numTraitsPL->get(PT); i++){
      genomeNames.insert("G:" + std::to_string(i) + "::");
      popFileColumns.push_back(std::to_string(i) + "::score");
    }

    
    // columns to be added to ave file
    popFileColumns.push_back("score");
    popFileColumns.push_back("score_VAR"); // specifies to also record the
                                         // variance (performed automatically
                                         // because _VAR)
}

void NKLandscapeWorld::evaluateSolo(std::shared_ptr<Organism> org, int analyze,int visualize, int debug) {
    auto N = N_PL->get(PT);
    auto K = K_PL->get(PT);

    double score = 0.0;

    for (auto &name : genomeNames){
        auto name2 = name.substr(2);
        auto genomeHandler = org->genomes[name2]->newHandler(org->genomes[name2], true);
        std::vector<int> siteBuffer = std::vector<int>();

        for (int i = 0; i < N; i++){
            auto value = genomeHandler->readInt(0, 1);
            siteBuffer.push_back(value);
    	}

        double genomeScore = 0;
        for (int i=0; i < N ; i++){
            auto index = 0;
            for (int j=0; j < K ; j++){
    	    	index += std::pow(2,j) * siteBuffer[(i+j)%N];
            }
            genomeScore += landscape[i][index]/N;
        }
	score += genomeScore;
        org->dataMap.append(name2 + "score",genomeScore);
    }

    org->dataMap.append("score", score);
    if (visualize)
        std::cout << "organism with ID " << org->ID << " scored " << score<< std::endl;  
}

void NKLandscapeWorld::evaluate(std::map<std::string, std::shared_ptr<Group>> &groups,int analyze, int visualize, int debug) {
  int popSize = groups[groupNamePL->get(PT)]->population.size();
  for (int i = 0; i < popSize; i++) {
    evaluateSolo(groups[groupNamePL->get(PT)]->population[i], analyze,visualize, debug);
  }
}

std::unordered_map<std::string, std::unordered_set<std::string>>
NKLandscapeWorld::requiredGroups() {
  return { {groupNamePL->get(PT), genomeNames } };
}
