//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#include "SawGenesWorld.h"

std::shared_ptr<ParameterLink<std::string>> SawGenesWorld::modePL =
Parameters::register_parameter("WORLD_SawGenes-mode", (std::string)"bit",
    "mode is either bit or int.\n"
    "  bit: genome will be read as (numGenes x geneSize) bits. Each set of geneSize bits, if bits are all within 1 of eachother, gene will converted using geneToScore, else gene score will be -100 (lethal). Agent score will be the sum of the gene scores.\n"
    "  int: genome will be read is numGenes ints each converted to a gene score using geneToScore. geneToScore will be used as a signal step in a repeating fitness function. genome alphabet size will be used to determin site value ranges. Agent score will be the sum of the gene scores.");

std::shared_ptr<ParameterLink<int>> SawGenesWorld::numGenesPL =
Parameters::register_parameter("WORLD_SawGenes-numGenes", 100,
    "how many genes?");


std::shared_ptr<ParameterLink<int>> SawGenesWorld::geneSizePL =
Parameters::register_parameter("WORLD_SawGenes-geneSize", 5,
    "how many sites in each gene?");

std::shared_ptr<ParameterLink<std::string>> SawGenesWorld::geneToScorePL =
Parameters::register_parameter("WORLD_SawGenes-geneToScore", (std::string)"0,-.1,-.2,-.3,-.4,10",
    "how to convert each gene to a gene score. agent score = sum of gene scores");

std::shared_ptr<ParameterLink<std::string>> SawGenesWorld::popSizePL =
Parameters::register_parameter("WORLD_SawGenes-popSize", (std::string)"-1x-1",
    "',' seperated list of 'x' seperated values. [popSize:generation,popSize:generation,...], times must be in order\n"
    "  e.g. 50x0,100x1000 will set popSize to 50 at generation 0, and then to 100 at 1000"
    "  use -1x-1 to use defualt popSize for all time"
    "  NOTE: the population genomes will be reset, even if the pop size does not change");


SawGenesWorld::SawGenesWorld(std::shared_ptr<ParametersTable> PT) : AbstractWorld(PT) {
    
    std::cout << "\nSetting up SawTooth World..." << std::endl;

    numGenes = numGenesPL->get(PT);
    geneSize = geneSizePL->get(PT);
    convertCSVListToVector(geneToScorePL->get(PT), geneToScore);

    std::cout << "    numGenes = " << numGenes << std::endl;
    std::cout << "    geneSize = " << geneSize << std::endl;
    std::cout << "    geneToScore = " << geneToScorePL->get(PT) << std::endl;

    if (geneToScore.size() != geneSize + 1) { // don't forget the value for 0
        std::cout << "      in SawGenesWorld constructor :: geneToScore size must == (geneSize + 1)... but it does not!" << std::endl;
        exit(1);
    }

    mode = modePL->get(PT);

    if (mode == "bit") {
        std::cout << "    mode = " << mode << std::endl;
        std::cout << "      genome should have " << numGenes * geneSize << " sites." << std::endl;
        std::cout << "      but bit mode is grabo... think think think... think. exiting..." << std::endl;
        exit(1);
    }
    else if (mode == "int") {
        std::cout << "    mode = " << mode << std::endl;
        std::cout << "      genome should have " << numGenes << " sites." << std::endl;
    }
    else {
        std::cout << "    unknown mode '" << mode << "'. exiting..." << std::endl;
        exit(1);
    }

    std::vector<std::string> temp1;
    std::vector<int> temp2;
    convertCSVListToVector(popSizePL->get(PT), temp1);

    if (temp1[0] != "-1x-1"){
        popSizes = {};
        popTimes = {};
        std::cout << "    found the following popSize instructions:" << std::endl;
        for (auto elem : temp1) {
            convertCSVListToVector(elem, temp2, 'x');
            popSizes.push_back(temp2[0]);
            popTimes.push_back(temp2[1]);
            std::cout << "      @ time " << popTimes.back() << " set population size to " << popSizes.back() << std::endl;
        }
        popSizes.push_back(-1);
        popTimes.push_back(-1);
    }

    // popFileColumns tell MABE what data should be saved to pop.csv files
	popFileColumns.clear();
    popFileColumns.push_back("score");
    popFileColumns.push_back("perfectGenesCount");
}

// the evaluate function gets called every generation. evaluate should set values on organisms datamaps
// that will be used by other parts of MABE for things like reproduction and archiving
auto SawGenesWorld::evaluate(std::map<std::string, std::shared_ptr<Group>>& groups, int analyze, int visualize, int debug) -> void {
    bool verbose = 0;

    int popSize = groups[groupName]->population.size();

    //std::cout << Global::update << " popTimes[" << popSizeIndex << "] = " << popTimes[popSizeIndex] << std::endl;
    //std::cout << Global::update << " popSizes[" << popSizeIndex << "] = " << popSizes[popSizeIndex] << std::endl;

    if (Global::update == 0 || Global::update == popTimes[popSizeIndex]) {
        // resize the population, if smaller then we can just kill the tail and resize
        //   if bigger we need to birth from last
        if (Global::update == popTimes[popSizeIndex] && popSizes[popSizeIndex] < popSize) { // shrink population
            //std::cout << "shirnk" << std::endl;
            for (int i = popSizes[popSizeIndex]; i < popSize; i++) {

                //std::cout << i << " " << popSize << " " << popSizes[popSizeIndex] << std::endl;

                groups[groupName]->population[i]->kill();
            }
            popSize = popSizes[popSizeIndex];
            groups[groupName]->population.resize(popSize);
        }
        else if (Global::update == popTimes[popSizeIndex] && popSizes[popSizeIndex] > popSize) { // grow population
            //std::cout << "growing" << std::endl;
            for (int i = popSize; i < popSizes[popSizeIndex]; i++) {
                
                //std::cout << i << " " << popSize << " " << popSizes[popSizeIndex] << std::endl;

                groups[groupName]->population.push_back(
                    groups[groupName]->population[0]->makeMutatedOffspringFrom(groups[groupName]->population[0])); // don't care who, genome is about to be reset
            }
            popSize = popSizes[popSizeIndex];
        }
        // else, popSize did not change, do nothing

        if (Global::update == popTimes[popSizeIndex]) {
            popSizeIndex++;
        }

        std::cout << "in SawTooth World :: initializing population... population size now: " << popSize << ", next update @ time " << popTimes[popSizeIndex] << " to " << popSizes[popSizeIndex] << std::endl;

        int sitesPerGene = (mode == "bit") ? geneSize : 1;
        for (auto& org : groups[groupName]->population) {
            auto gh = org->genomes[genomeName]->newHandler(org->genomes[genomeName]);
            for (int i = 0; i < numGenes; i++){
                for (int j = 0; j < sitesPerGene; j++) {
                    gh->writeInt((int)(org->genomes[genomeName]->getAlphabetSize() / 2), 0, org->genomes[genomeName]->getAlphabetSize() - 1);
                }
            }
        }
    }

    // in this world, organisms donot interact, so we can just iterate over the population
    for (int i = 0; i < popSize; i++) {

        // create a shortcut to access the organism and organisms brain
        auto org = groups[groupName]->population[i];
        auto genome = org->genomes[genomeName];

        if (0) { // if you want to see the bit version of this genome
            auto gh = genome->newHandler(genome);
            for (int i = 0; i < numGenes; i++) {
                for (int j = 0; j < geneSize; j++) {
                    std::cout << gh->readInt(0, 1);
                }
                std::cout << " ";
            }
            std::cout << std::endl;
        }

        auto genomeHandler = genome->newHandler(genome);

        double score = 0;
        std::vector<int> geneValues(geneSize);
        int minGeneValue;
        int perfectGenesCount = 0;
        bool badGene = false;

        for (int i = 0; i < numGenes; i++) {
            if (mode == "bit") {
                int geneSum = 0;
                for (int j = 0; j < geneSize; j++) { // collect values for this gene and sum of values
                    geneValues[j] = genomeHandler->readInt(0, genome->getAlphabetSize() - 1);
                    geneSum += geneValues[j];
                }
                int minGeneValue = *std::min_element(geneValues.begin(), geneValues.end()); // smallest site value in gene
                int maxGeneValue = *std::max_element(geneValues.begin(), geneValues.end()); // greatest site vlaue in gene
                if (maxGeneValue - minGeneValue > 1) {
                    badGene = true;
                    //this gene does not add to the score
                    if (verbose) {
                        std::cout << "leathal gene" << std::endl;
                    }
                }
                else {
                    score += geneToScore[geneSize] * minGeneValue + geneToScore[std::accumulate(geneValues.begin(), geneValues.end(), -1 * minGeneValue * geneSize)];
                    if (verbose) {
                        std::cout << score << "  " << minGeneValue << " : " << geneValues[0] << "," << geneValues[1] << "," << geneValues[2] << "," << geneValues[3] << "," <<
                            geneValues[4] << "," << geneValues[5] << "," << geneValues[6] << "," << geneValues[7] << " -> " <<
                            std::accumulate(geneValues.begin(), geneValues.end(), -1 * minGeneValue * geneSize) << " = " <<
                            geneToScore[std::accumulate(geneValues.begin(), geneValues.end(), -1 * minGeneValue * geneSize)] << std::endl;
                    }
                }
                perfectGenesCount += minGeneValue;
            }
            else { // int mode
                double siteValue = genomeHandler->readInt(0, genome->getAlphabetSize() - 1);
                score += geneToScore[geneSize] * (int)(siteValue / geneSize) + geneToScore[loopMod(siteValue, geneSize)];
                perfectGenesCount += (int)(siteValue / geneSize);
                if (verbose) {
                    std::cout << siteValue << "  " << score << "  " << perfectGenesCount << std::endl;
                }
            }
        }
        if (badGene) {
            score = 0;
        }
        org->dataMap.append("score", score);
        org->dataMap.append("perfectGenesCount", perfectGenesCount);
        if (verbose) {
            std::cout << "   score: " << score << "   perfectGenesCount: " << perfectGenesCount << std::endl;
        }
    }
}

// the requiredGroups function lets MABE know how to set up populations of organisms that this world needs
auto SawGenesWorld::requiredGroups() -> std::unordered_map<std::string, std::unordered_set<std::string>> {
	return { { groupName, { "G:"+genomeName } } };
}
