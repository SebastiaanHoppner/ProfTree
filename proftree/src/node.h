#include <Rcpp.h>
#include <iostream>
#include <R.h>
#include <Rinternals.h>
#include "variable.h"
#include "empChurnCpp.h"
#define pi 3.14159

using namespace std;

class Node{
public:
  int pos;
  int* splitV;
  double* splitP;
  int** csplit;
  Node* leftChild;
  Node* rightChild;
  int* nInst;
  int* nVar;
  int* localClassification;
  double** data;
  int* nClassesDependendVar;
  variable **variables;
  
  double performanceLeftTerminal;
  double performanceRightTerminal;
  
  int sumLocalWeights; // weighted number of instances
  int sumLeftLocalWeights;
  int sumRightLocalWeights;
  
  double predictionInternalNode;
  double predictionLeftTerminal;
  double predictionRightTerminal;
  
  Node(int splitN, int* splitV, double* splitP, int** csplit, Node* leaftChild, Node* rightChild, double** data, int* nInst, int* nVar, variable** variables);
  ~Node();
  int partition( int* classification, int* weights, variable** variables, int* nNodes, int minbucket, int minsplit);
  void calculateChildNodeScore(bool leftNode, int* weights, IntegerVector& nObsClass1PerNode, IntegerVector& nObsPerNode);
};
