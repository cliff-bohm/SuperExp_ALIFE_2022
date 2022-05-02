//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#include "DriftPoolRouletteOptimizer.h"

std::shared_ptr<ParameterLink<int>> DriftPoolRouletteOptimizer::numberParentsPL =
	Parameters::register_parameter("OPTIMIZER_DRIFTPOOLROULETTE-numberParents", 1, "number of parents used to produce offspring");

std::shared_ptr<ParameterLink<int>> DriftPoolRouletteOptimizer::driftPoolSizePL =
	Parameters::register_parameter("OPTIMIZER_DRIFTPOOLROULETTE-driftPoolSize", 2, "Number of organisms in the drift pool.");

std::shared_ptr<ParameterLink<double>> DriftPoolRouletteOptimizer::poolDecayRatePL =
	Parameters::register_parameter("OPTIMIZER_DRIFTPOOLROULETTE-poolDecayRate", 0.01, "Per organism per generation probability to decay and be replaced by max");


std::shared_ptr<ParameterLink<std::string>> DriftPoolRouletteOptimizer::optimizeValuePL =
Parameters::register_parameter("OPTIMIZER_DRIFTPOOLROULETTE-optimizeValue", (std::string) "DM_AVE[score]", "value to optimize (MTree)");

std::shared_ptr<ParameterLink<std::string>> DriftPoolRouletteOptimizer::remapFunctionPL =
Parameters::register_parameter("OPTIMIZER_DRIFTPOOLROULETTE-remapFunction", (std::string) "POW[2.72,$optVal$-$maxOptVal$]", "remap optimizeValue to affect strength of selection\n"
	"uses MTree, but adds the following options that can be used in place of MTree functions:\n"
	"$optVal$ = score of current organism\n"
	"$maxOptVal$ = maximum score in population\n"
	"$minOptVal$ = minimum score in population\n"
    "$aveOptVal$ = average score in population\n"
	"if NONE, no remap is preformed.");

std::shared_ptr<ParameterLink<std::string>> DriftPoolRouletteOptimizer::selectionMethodPL =
Parameters::register_parameter("OPTIMIZER_DRIFTPOOLROULETTE-selectionMethod", (std::string) "roulette", "selectionMethod roulette or tournament:[size]");

// --------------------------------------------------------------------------------------------
//CONSTRUCTOR ==============================================================================
// --------------------------------------------------------------------------------------------
DriftPoolRouletteOptimizer::DriftPoolRouletteOptimizer(std::shared_ptr<ParametersTable> PT_)
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
		std::cout << "  in DriftPoolRouletteOptimizer: found a bad selection method '" << selectionMethodPL->get(PT) << "'" << std::endl;
	}
}

// --------------------------------------------------------------------------------------------
// HELPER FUNCTIONS ============================================================================
// --------------------------------------------------------------------------------------------
int DriftPoolRouletteOptimizer::selectParentRoulette(std::vector<double>& scores, double maxScore, double minScore, int popSize) {

	if (maxScore <= 0 || maxScore == minScore) {
		return Random::getIndex(popSize);
	}
	int parent;
	do {
		parent = Random::getIndex(popSize);
	} while (!Random::P(scores[parent] / maxScore));
	return parent;
}

int DriftPoolRouletteOptimizer::selectParentTournament(std::vector<double>& scores, double maxScore, double minScore, int popSize, double tournamentSize) {
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

void DriftPoolRouletteOptimizer::stringReplace(std::string& s, const std::string& search, const std::string& replace) {
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
void DriftPoolRouletteOptimizer::optimize(std::vector<std::shared_ptr<Organism>>& population) {
    	std::cout << std::endl;
	auto popSize = population.size();

	std::vector<double> scores(popSize, 0);
	std::vector<double> remappedScores(popSize, 0);
	double aveScore = 0;
	double maxScore = optimizeValueMT->eval(population[0]->dataMap, PT)[0];
	double minScore = maxScore;

	killList.clear();

	for (size_t i = 0; i < popSize; i++) {
		killList.insert(population[i]);
		double opVal = optimizeValueMT->eval(population[i]->dataMap, PT)[0];
		scores[i] = opVal;
		aveScore += opVal;
		population[i]->dataMap.set("optimizeValue", opVal);
        maxScore = std::max(maxScore, opVal);
		minScore = std::min(minScore, opVal);
	}

	aveScore /= popSize;
    
	std::vector<double> maxPoolRemappedScores(1);
	std::vector<std::vector<double>> remapVect({ { minScore, aveScore, maxScore, 0 } });
	if (doRemap) {
		for (size_t i = 0; i < popSize; i++) {
			remapVect[0][3] = scores[i];
			remappedScores[i] = remapFunctionMT->eval(population[i]->dataMap, PT, remapVect)[0];
			population[i]->dataMap.set("remappedOptimizeValue", remappedScores[i]);
		}
		remapVect[0][3] = maxPoolScores[0];
		maxPoolRemappedScores[0] = remapFunctionMT->eval(maxPool[0]->dataMap, PT, remapVect)[0];;
	}
	else {
		remappedScores = scores;
		maxPoolRemappedScores[0] = maxPoolScores[0];
	}

	int remappedScoresMaxIndex = std::max_element(remappedScores.begin(), remappedScores.end())-remappedScores.begin();
	//int remappedScoresMaxIndex_poponly = std::max_element(remappedScores.begin(), remappedScores.end()-driftPoolSize)-remappedScores.begin();

	double remappedScoresMax = remappedScores[remappedScoresMaxIndex];
	//double remappedScoresMax_poponly = remappedScores[remappedScoresMaxIndex_poponly];

	double remappedScoresMin = *std::min_element(remappedScores.begin(), remappedScores.end());
	//double remappedScoresMin_poponly = *std::min_element(remappedScores.begin(), remappedScores.end()-driftPoolSize);


	if (Global::update == 0){
		for (auto o : population) {
			o->dataMap.set("pool", -1);
		}
	}
	if (Global::update == 0 || (remappedScoresMax > maxPoolRemappedScores[0])) {
		maxPoolScores[0] = scores[remappedScoresMaxIndex];
		if (Global::update > 0) { // if 0 then maxPool[0] will be nil
			maxPool[0]->kill();
			//std::cout << "parent count: " << maxPool[0]->parents[0]->offspringCount << std::endl;
		}
		//maxPool[0] = population[remappedScoresMaxIndex]->makeCopy(population[remappedScoresMaxIndex]->PT);
		maxPool[0] = std::make_shared<Organism>(population[remappedScoresMaxIndex]->PT);
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
			population.back()->dataMap.set("pool", 1);
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
			population.back()->dataMap.set("pool", 1);
		}
	}

	int replacementCount = 0;
	for (int i = popSize - driftPoolSize; i < popSize ; i++){
	    if (Random::P(decayRate)){
	        population.push_back(maxPool[0]->makeMutatedOffspringFrom(maxPool[0]));
			//std::cout << "\t\ta pool org (" << i << ") was reset to max( " << remappedScoresMaxIndex_poponly << " )." << std::endl;
			population.back()->dataMap.set("pool", 0);
			replacementCount++;
		}
	    else {
	        population.push_back(population[i]->makeMutatedOffspringFrom(population[i]));
			population.back()->dataMap.set("pool", 0);
		}
	}

	for (int i = 0; i < population.size(); i++) {
		if (population[i]->timeOfBirth < Global::update) {
			population[i]->dataMap.set("roulette_numOffspring", population[i]->offspringCount);
		}
	}
	maxPool[0]->dataMap.set("roulette_numOffspring", maxPool[0]->offspringCount);


	std::cout << "maxPool score = " << std::to_string(maxPoolScores[0]) << "   " << replacementCount << " replacements." << std::endl;
	std::cout << "   max = " << std::to_string(maxScore) << "   ave = " << std::to_string(aveScore) << "   min = " << std::to_string(minScore);
}
