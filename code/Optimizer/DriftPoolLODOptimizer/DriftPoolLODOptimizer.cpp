//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#include "DriftPoolLODOptimizer.h"

std::shared_ptr<ParameterLink<int>> DriftPoolLODOptimizer::numberParentsPL =
	Parameters::register_parameter("OPTIMIZER_DRIFTPOOLLOD-numberParents", 1, "number of parents used to produce offspring");

std::shared_ptr<ParameterLink<int>> DriftPoolLODOptimizer::driftPoolSizePL =
	Parameters::register_parameter("OPTIMIZER_DRIFTPOOLLOD-driftPoolSize", 2, "Number of organisms in the drift pool.");

std::shared_ptr<ParameterLink<double>> DriftPoolLODOptimizer::poolDecayRatePL =
	Parameters::register_parameter("OPTIMIZER_DRIFTPOOLLOD-poolDecayRate", 0.01, "Per organism per generation probability to decay and be replaced by max");


std::shared_ptr<ParameterLink<std::string>> DriftPoolLODOptimizer::optimizeValuePL =
Parameters::register_parameter("OPTIMIZER_DRIFTPOOLLOD-optimizeValue", (std::string) "DM_AVE[score]", "value to optimize (MTree)");

std::shared_ptr<ParameterLink<std::string>> DriftPoolLODOptimizer::remapFunctionPL =
Parameters::register_parameter("OPTIMIZER_DRIFTPOOLLOD-remapFunction", (std::string) "POW[2.72,$optVal$-$maxOptVal$]", "remap optimizeValue to affect strength of selection\n"
	"uses MTree, but adds the following options that can be used in place of MTree functions:\n"
	"$optVal$ = score of current organism\n"
	"$maxOptVal$ = maximum score in population\n"
	"$minOptVal$ = minimum score in population\n"
    "$aveOptVal$ = average score in population\n"
	"if NONE, no remap is preformed.");

std::shared_ptr<ParameterLink<std::string>> DriftPoolLODOptimizer::selectionMethodPL =
Parameters::register_parameter("OPTIMIZER_DRIFTPOOLLOD-selectionMethod", (std::string) "roulette", "selectionMethod roulette or tournament:[size]");

std::shared_ptr<ParameterLink<std::string>> DriftPoolLODOptimizer::decayModePL =
Parameters::register_parameter("OPTIMIZER_DRIFTPOOLLOD-decayMode", (std::string) "GlobalMax", "how to regenerate decayed Super Explorers: GlobalMax,PopMax,PopProp");

// --------------------------------------------------------------------------------------------
//CONSTRUCTOR ==============================================================================
// --------------------------------------------------------------------------------------------
DriftPoolLODOptimizer::DriftPoolLODOptimizer(std::shared_ptr<ParametersTable> PT_)
	: AbstractOptimizer(PT_) {

	numberParents = numberParentsPL->get(PT);

	optimizeValueMT = stringToMTree(optimizeValuePL->get(PT));
    
    driftPoolSize = driftPoolSizePL->get(PT);
    decayRate = poolDecayRatePL->get(PT);
    maxPool.resize(1);
    maxPoolScores.resize(1);

	if (remapFunctionPL->get(PT) == "NONE") {
		doRemap = false;
	}
	else {
		doRemap = true;
		std::string remapString = remapFunctionPL->get(PT);
		stringReplace(remapString, "$minOptVal$", "VECT[0,0]");
		stringReplace(remapString, "$aveOptVal$", "VECT[0,1]");
		stringReplace(remapString, "$maxOptVal$", "VECT[0,2]");
		stringReplace(remapString, "$optVal$", "VECT[0,3]");
		remapFunctionMT = stringToMTree(remapString);
	}
	popFileColumns.clear();
	popFileColumns.push_back("optimizeValue");
	//if (doRemap) {
	//	popFileColumns.push_back("remappedOptimizeValue");
	//}
    if (doRemap) {
        popFileColumns.push_back("remappedOptimizeValue");
        optimizeFormula = stringToMTree("DM_AVE[remappedOptimizeValue]");
    }
    else {
        optimizeFormula = stringToMTree("DM_AVE[optimizeValue]");
    }

	std::vector<std::string> tempList;
	convertCSVListToVector(selectionMethodPL->get(PT), tempList, ':');
	if (tempList[0] == "roulette") {
		selectionMethod = "roulette";
	}
	else if (tempList[0] == "tournament") {
		selectionMethod = "tournament";
		convertString(tempList[1], tournamentSize);
	}
	else {
		std::cout << "  in DriftPoolLODOptimizer: found a bad selection method '" << selectionMethodPL->get(PT) << "'" << std::endl;
	}

	if (decayModePL->get(PT) == "GlobalMax"){
		std::cout << "Decay Mode: GlobalMax"<< std::endl;
	}
	else if (decayModePL->get(PT) == "PopMax"){
		std::cout << "Decay Mode: PopMax"<< std::endl;
	}
	else if (decayModePL->get(PT) == "PopProp"){
		std::cout << "Decay Mode: PopProp"<< std::endl;
	}
	else{
		std::cout << "  in DrfotPoolLODOptimizer: found a bad decay mode '" << decayModePL->get(PT) << "'" << std::endl;
	}
}

// --------------------------------------------------------------------------------------------
// HELPER FUNCTIONS ============================================================================
// --------------------------------------------------------------------------------------------
int DriftPoolLODOptimizer::selectParentRoulette(std::vector<double>& scores, double maxScore, double minScore, int popSize) {

	if (maxScore <= 0 || maxScore == minScore) {
		return Random::getIndex(popSize);
	}
	int parent;
	do {
		parent = Random::getIndex(popSize);
	} while (!Random::P(scores[parent] / maxScore));
	return parent;
}

int DriftPoolLODOptimizer::selectParentTournament(std::vector<double>& scores, double maxScore, double minScore, int popSize, double tournamentSize) {
	if (maxScore <= 0 || maxScore == minScore) {
		return Random::getIndex(popSize);
	}
	int actualTournamentSize = int(tournamentSize);
	if (tournamentSize != actualTournamentSize) {
		actualTournamentSize += (Random::P(tournamentSize - actualTournamentSize)) ? 1 : 0;
	}
	int parent = Random::getIndex(popSize);
	int challenger;
	for (int i = 1; i < actualTournamentSize; i++){
		challenger = Random::getIndex(popSize);
		if (scores[challenger] > scores[parent]) {
			parent = challenger;
		}
	}
	return parent;
}

void DriftPoolLODOptimizer::stringReplace(std::string& s, const std::string& search, const std::string& replace) {
	for (size_t pos = 0; ; pos += replace.length()) {
		// Locate the substd::string to replace
		pos = s.find(search, pos);
		if (pos == std::string::npos) break;
		// Replace by erasing and inserting
		s.erase(pos, search.length());
		s.insert(pos, replace);
	}
}

// --------------------------------------------------------------------------------------------
//OPTIMIZE ====================================================================================
// --------------------------------------------------------------------------------------------
void DriftPoolLODOptimizer::optimize(std::vector<std::shared_ptr<Organism>>& population) {
    	std::cout << std::endl; //output formatting

	//local variables
	auto popSize = population.size();
	std::vector<double> scores(popSize, 0);
	std::vector<double> remappedScores(popSize, 0);
	double aveScore = 0;
	double maxScore = optimizeValueMT->eval(population[0]->dataMap, PT)[0];
	double minScore = maxScore;

	//clear kill list
	killList.clear();
	
	//for every organism...
	for (size_t i = 0; i < popSize; i++) {
		// add to kill list (all die by default)
		killList.insert(population[i]);

		//retreive optimize value from datamap
		double opVal = optimizeValueMT->eval(population[i]->dataMap, PT)[0];
		
		// collect scores into a vector, update ave, min, and max.
		scores[i] = opVal;
		aveScore += opVal;
		population[i]->dataMap.set("optimizeValue", opVal);//set optimizeValue
        	maxScore = std::max(maxScore, opVal);
		minScore = std::min(minScore, opVal);
	}
	// complete average calculation
	aveScore /= popSize;
    
	//more local variables... required min, ave, and max.
	std::vector<double> maxPoolRemappedScores(1);
	std::vector<std::vector<double>> remapVect({ { minScore, aveScore, maxScore, 0 } }); //this acts like a data packet. Slot 3 is reused.

	//if there's a remap function, run everything through the MTree.
	if (doRemap) {
		for (size_t i = 0; i < popSize; i++) {
			remapVect[0][3] = scores[i]; //the score to be remapped is in slot 3
			remappedScores[i] = remapFunctionMT->eval(population[i]->dataMap, PT, remapVect)[0];
			population[i]->dataMap.set("remappedOptimizeValue", remappedScores[i]);
		}
		remapVect[0][3] = maxPoolScores[0]; //we need to remap max too in case it would be diferent this time.
		maxPoolRemappedScores[0] = remapFunctionMT->eval(maxPool[0]->dataMap, PT, remapVect)[0];;
	}
	else {
		remappedScores = scores;
		maxPoolRemappedScores[0] = maxPoolScores[0];
	}

	// get the (potentially) new min and max after the remap.
	int remappedScoresMaxIndex = std::max_element(remappedScores.begin(), remappedScores.end())-remappedScores.begin();
	double remappedScoresMax = remappedScores[remappedScoresMaxIndex];
	double remappedScoresMin = *std::min_element(remappedScores.begin(), remappedScores.end());
	
	//int remappedScoresMaxIndex_poponly = std::max_element(remappedScores.begin(), remappedScores.end()-driftPoolSize)-remappedScores.begin();
	//double remappedScoresMax_poponly = remappedScores[remappedScoresMaxIndex_poponly];
	//double remappedScoresMin_poponly = *std::min_element(remappedScores.begin(), remappedScores.end()-driftPoolSize);


	// if this is the first update, mark the source as -1 for honest accounting of LOD.
	// also store the first max in the max pool
	if (Global::update == 0){
		for (auto o : population) {
			o->dataMap.set("pool", -1);
		}
		// actually set the max pool's score to the new max
		maxPoolScores[0] = scores[remappedScoresMaxIndex];
		
		//store the organism in the max pool
		//requires a manual copy so LOD is preserved.
		maxPool[0] = std::make_shared<Organism>(population[remappedScoresMaxIndex]->PT); //make a NEW organism based on the PT; NOT a shared pointer to organism
		for (auto genome : population[remappedScoresMaxIndex]->genomes) {
			maxPool[0]->genomes[genome.first] = genome.second->makeCopy(genome.second->PT);
		}
		for (auto brain : population[remappedScoresMaxIndex]->brains) {
			maxPool[0]->brains[brain.first] = brain.second->makeCopy(brain.second->PT);
		}
		maxPool[0]->dataMap = population[remappedScoresMaxIndex]->dataMap;

		maxPool[0]->offspringCount = 0;
		
		maxPool[0]->parents = { population[remappedScoresMaxIndex] };
		population[remappedScoresMaxIndex]->offspringCount++;
		
		maxPool[0]->ancestors = population[remappedScoresMaxIndex]->ancestors;
		maxPool[0]->timeOfBirth = Global::update;
		maxPool[0]->timeOfDeath = -1;
		maxPool[0]->alive = 1;

		maxPool[0]->dataMap.set("pool", 2);
		maxPool[0]->trackOrganism = true;
	}

	
	// Do selection. Organisms are chosen from both the population and drift pool. (same vector, diferent indices)
	// Pop(t-1),DP(t-1) -> pop(t)
	std::shared_ptr<Organism> parent; // used for single parent
	std::vector<std::shared_ptr<Organism>> parents; // used for multi-parent

	for (int i = 0; i < popSize-driftPoolSize; i++) {
		if (numberParents == 1) {
			if (selectionMethod == "roulette") {
				parent = population[selectParentRoulette(remappedScores, remappedScoresMax, remappedScoresMin, popSize)];
			}
			else {
				parent = population[selectParentTournament(remappedScores, remappedScoresMax, remappedScoresMin, popSize, tournamentSize)];
			}
			population.push_back(parent->makeMutatedOffspringFrom(parent)); // add to population
		}
		else { // numberParents != 1
			parents.clear();
			do {
				if (selectionMethod == "roulette") {
					parents.push_back(population[selectParentRoulette(remappedScores, remappedScoresMax, remappedScoresMin, popSize)]);
				}
				else {
					parents.push_back(population[selectParentTournament(remappedScores, remappedScoresMax, remappedScoresMin, popSize, tournamentSize)]);
				}
			} while (static_cast<int>(parents.size()) < numberParents);
			population.push_back(parents[0]->makeMutatedOffspringFromMany(parents)); // push to population
		}
		population.back()->dataMap.set("pool", 1);
	}
	

	//decay and replace any drift pool organisms according to the parameters.
	int replacementCount = 0;
	for (int i = popSize - driftPoolSize; i < popSize ; i++){
	    if (Random::P(decayRate)){
			if (decayModePL->get(PT) == "GlobalMax"){
				population.push_back(maxPool[0]->makeMutatedOffspringFrom(maxPool[0]));
			}
			else if (decayModePL->get(PT) == "PopMax"){
				population.push_back(population[remappedScoresMaxIndex]->makeMutatedOffspringFrom(population[remappedScoresMaxIndex]));
			}
			else if (decayModePL->get(PT) == "PopProp"){
				if (numberParents == 1) {
					if (selectionMethod == "roulette") {
						parent = population[selectParentRoulette(remappedScores, remappedScoresMax, remappedScoresMin, popSize)];
					}
					else {
						parent = population[selectParentTournament(remappedScores, remappedScoresMax, remappedScoresMin, popSize, tournamentSize)];
					}
					population.push_back(parent->makeMutatedOffspringFrom(parent)); // add to population
				}
				else { // numberParents != 1
					parents.clear();
					do {
						if (selectionMethod == "roulette") {
							parents.push_back(population[selectParentRoulette(remappedScores, remappedScoresMax, remappedScoresMin, popSize)]);
						}
						else {
							parents.push_back(population[selectParentTournament(remappedScores, remappedScoresMax, remappedScoresMin, popSize, tournamentSize)]);
						}
					} while (static_cast<int>(parents.size()) < numberParents);
					population.push_back(parents[0]->makeMutatedOffspringFromMany(parents)); // push to population
				}
			}

			population.back()->dataMap.set("pool", 0);
			replacementCount++;
		}
		else {
	        	population.push_back(population[i]->makeMutatedOffspringFrom(population[i]));
			population.back()->dataMap.set("pool", 0);
		}
	}
	

	//update offspring counts
	for (int i = 0; i < population.size(); i++) {
		if (population[i]->timeOfBirth < Global::update) {
			population[i]->dataMap.set("roulette_numOffspring", population[i]->offspringCount);
		}
	}
	maxPool[0]->dataMap.set("roulette_numOffspring", maxPool[0]->offspringCount);

	
	//print to console
	std::cout << "maxPool score = " << std::to_string(maxPoolScores[0]) << "   " << replacementCount << " replacements." << std::endl;
	std::cout << "   max = " << std::to_string(maxScore) << "   ave = " << std::to_string(aveScore) << "   min = " << std::to_string(minScore);
	

	// any time max pool needs to be updated...
	if (remappedScoresMax > maxPoolRemappedScores[0]) {
		// actually set the max pool's score to the new max
		maxPoolScores[0] = scores[remappedScoresMaxIndex];
		
		//kill the old max so it can be free	
		maxPool[0]->kill();

		//store the organism in the max pool
		//requires a manual copy so LOD is preserved.
		maxPool[0] = std::make_shared<Organism>(population[remappedScoresMaxIndex]->PT); //make a NEW organism based on the PT; NOT a shared pointer to organism
		for (auto genome : population[remappedScoresMaxIndex]->genomes) {
			maxPool[0]->genomes[genome.first] = genome.second->makeCopy(genome.second->PT);
		}
		for (auto brain : population[remappedScoresMaxIndex]->brains) {
			maxPool[0]->brains[brain.first] = brain.second->makeCopy(brain.second->PT);
		}
		maxPool[0]->dataMap = population[remappedScoresMaxIndex]->dataMap;

		maxPool[0]->offspringCount = 0;
		
		maxPool[0]->parents = { population[remappedScoresMaxIndex] };
		population[remappedScoresMaxIndex]->offspringCount++;
		
		maxPool[0]->ancestors = population[remappedScoresMaxIndex]->ancestors;
		maxPool[0]->timeOfBirth = Global::update;
		maxPool[0]->timeOfDeath = -1;
		maxPool[0]->alive = 1;

		maxPool[0]->dataMap.set("pool", 2);
		maxPool[0]->trackOrganism = true;
	}
}
