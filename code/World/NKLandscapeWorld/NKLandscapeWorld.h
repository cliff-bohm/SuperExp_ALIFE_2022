//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#pragma once

#include <World/AbstractWorld.h>

#include <cstdlib>
#include <thread>
#include <vector>

class NKLandscapeWorld : public AbstractWorld {

public:
  static std::shared_ptr<ParameterLink<int>> modePL;
  static std::shared_ptr<ParameterLink<int>> N_PL;
  static std::shared_ptr<ParameterLink<int>> K_PL;
  static std::shared_ptr<ParameterLink<int>> numTraitsPL;

  static std::shared_ptr<ParameterLink<std::string>> groupNamePL;
  static std::shared_ptr<ParameterLink<std::string>> genomeNamePL;

  NKLandscapeWorld(std::shared_ptr<ParametersTable> PT_ = nullptr);
  virtual ~NKLandscapeWorld() = default;

  std::vector<std::vector<double>> landscape;
  std::unordered_set<std::string> genomeNames;
  
  void evaluateSolo(std::shared_ptr<Organism> org, int analyze,int visualize, int debug);
  void evaluate(std::map<std::string, std::shared_ptr<Group>> &groups,int analyze, int visualize, int debug);

  virtual std::unordered_map<std::string, std::unordered_set<std::string>> requiredGroups() override;
};

