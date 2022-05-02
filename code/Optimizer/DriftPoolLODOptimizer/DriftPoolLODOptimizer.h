//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#include "../AbstractOptimizer.h"
#include "../../Utilities/MTree.h"

#include <iostream>
#include <sstream>

class DriftPoolLODOptimizer : public AbstractOptimizer {
public:
	static std::shared_ptr<ParameterLink<int>> numberParentsPL;
	static std::shared_ptr<ParameterLink<std::string>> optimizeValuePL;
	static std::shared_ptr<ParameterLink<std::string>> remapFunctionPL;
        static std::shared_ptr<ParameterLink<int>> driftPoolSizePL;
	static std::shared_ptr<ParameterLink<double>> poolDecayRatePL;
	static std::shared_ptr<ParameterLink<std::string>> selectionMethodPL;
	static std::shared_ptr<ParameterLink<std::string>> decayModePL;
	std::string selectionMethod;
	double tournamentSize;

	int selectParentRoulette(std::vector<double>& scores, double maxScore, double minScore, int popSize);
	int selectParentTournament(std::vector<double>& scores, double maxScore, double minScore, int popSize, double tournamentSize);

	void stringReplace(std::string& s, const std::string& search, const std::string& replace);

	int numberParents;
	int driftPoolSize;
	double decayRate;
	std::shared_ptr<Abstract_MTree> optimizeValueMT;
	std::shared_ptr<Abstract_MTree> remapFunctionMT;
	bool doRemap;
    std::vector<std::shared_ptr<Organism>> maxPool;
    std::vector<double> maxPoolScores;

	DriftPoolLODOptimizer(std::shared_ptr<ParametersTable> PT_ = nullptr);

	virtual void optimize(std::vector<std::shared_ptr<Organism>> &population) override;
};
