#include "tree.h"
#include "empChurnCpp.h"

Tree::Tree(int* nInstances, int* nVariables, double** data, int* weights, int *splitV, double *splitP, int** csplit, int* maxCat, int* nNodes, variable** variables, int* maxNode){
  // constructor used by crossover
  this->nInstances = nInstances;
  this->nVariables = nVariables;
  this->nNodes = *nNodes;
  this->maxNode = maxNode;
  this->maxCat = maxCat;
  this->splitV = new int[*this->maxNode];
  this->splitP = new double[*this->maxNode];
  this->variables = variables;
  this->nodes = new Node*[*this->maxNode];
  this->classification = new int[*this->nInstances];
  this->data = data;
  this->performance = 999999;
  this->csplit = new int*[*this->maxCat];
  this->weights = weights;
  
  for (int i = 0; i < *this->maxCat; i++)
    this->csplit[i] = new int[(*this->maxNode)];
  for (int v = 0; v < *this->maxNode; v++){
    for(int i = 0; i < *this->maxCat; i++){
      this->csplit[i][v] = csplit[i][v];
    }
    this->splitV[v] = splitV[v];
    this->splitP[v] = splitP[v];
  }
  for (int nodeNumber = *this->maxNode-1; nodeNumber >= 0; nodeNumber--){
    this->nodes[nodeNumber] = NULL;
    this->initNode(nodeNumber);
  }
  
} // end Tree


Tree::Tree(int* nInstances, int* nVariables, double** data, int* weights, int* maxCat, variable** variables, int* maxNode, int* minBucket, int* minSplit){
  // initializes a tree with a random root node
  this->nInstances = nInstances;
  this->nVariables = nVariables;
  this->nNodes = 1;
  this->maxNode = maxNode;
  this->maxCat = maxCat;
  this->splitV = new int[*this->maxNode];
  this->splitP = new double[*this->maxNode];
  this->variables = variables;
  this->nodes = new Node*[*this->maxNode];
  this->classification = new int[*this->nInstances];
  this->data = data;
  this->performance = 999999;
  this->csplit = new int*[*this->maxCat];
  this->weights = weights;
  
  for (int i = 0; i < *this->maxCat; i++)
    this->csplit[i] = new int[(*this->maxNode)];
  for (int v = 0; v < *this->maxNode; v++){
    for(int i = 0; i < *this->maxCat; i++){
      this->csplit[i][v] = 2;
    }
    this->splitV[v] = -999999;
    this->splitP[v] = -999999;
    this->nodes[v] = NULL;
  }
  
  this->splitV[0] = Tree::getUnifRandNumber(*this->nVariables-1);
  this->nodes[0] = NULL;
  this->initNode(0);
  
  if(variables[this->splitV[0]]->isCat == false){
    if((this->variables[this->splitV[0]]->nCats-1) > 1)
      this->splitP[0] = variables[this->splitV[0]]->sortedValues[getUnifRandNumber(this->variables[this->splitV[0]]->nCats-1)+1];
    else
      this->splitP[0] = variables[this->splitV[0]]->sortedValues[0];
  }else{
    this->randomizeCategories(0);
  }
  
  for(int i = 0; i <= 5000 && this->predictClass(*minBucket, *minSplit, false, 0) == false; i++){
    this->splitV[0] = getUnifRandNumber(*this->nVariables-1);
    if(variables[this->splitV[0]]->isCat == false){
      if((this->variables[this->splitV[0]]->nCats-1) > 1 )
        this->splitP[0] = variables[this->splitV[0]]->sortedValues[getUnifRandNumber(this->variables[this->splitV[0]]->nCats-1)+1];
      else
        this->splitP[0] = variables[this->splitV[0]]->sortedValues[0];
    }else{
      this->randomizeCategories(0);
    }
    if(i == 5000)
      this->splitV[0] = -999999;
  }
} // end Tree


Tree::~Tree(){
  // deconstructor destroys the instance and frees up memory
  for (int i = 0; i < *this->maxNode; i++){
    delete nodes[i];
  }
  delete [] nodes;
  nodes = NULL;
  delete [] classification;
  classification = NULL;
  delete [] splitP;
  splitP = NULL;
  delete [] splitV;
  splitV = NULL;
  for (int i = 0; i < *this->maxCat; i++)
    delete [] csplit[i];
  delete [] csplit;
  csplit = NULL;
  variables = NULL;
  data = NULL;
  maxNode = NULL;
  maxCat = NULL;
  nInstances = NULL;
  nVariables = NULL;
  weights = NULL;
} // end ~Tree


int Tree::getUnifRandNumber(int numberDistinctValues){
  return ((int)floorf(unif_rand()*((double)numberDistinctValues)))%numberDistinctValues; // % for the case unif_rand gives exactly 1
}


void Tree::initNode(int nodeNumber){
  // initializes a node
  if(this->splitV[nodeNumber] >= 0 && nodeNumber >= 0 ){
    int leftChild = -1;
    int rightChild = -1;
    // is leaf node?
    if( nodeNumber*2+2 < *this->maxNode){
      if( (this->splitV[nodeNumber*2+1]  ) >=  0 ){
        leftChild = nodeNumber*2+1;
      }
      if( (this->splitV[nodeNumber*2+2]) >=  0){
        rightChild = nodeNumber*2+2;
      }
    }
    
    if( leftChild <= 0 && rightChild <= 0){
      this->nodes[nodeNumber] = new Node(nodeNumber, &this->splitV[nodeNumber], &this->splitP[nodeNumber], this->csplit, NULL, NULL,
                                         this->data, this->nInstances, this->nVariables,  this->variables);
    }else if( leftChild <= 0 ){
      this->nodes[nodeNumber] = new Node(nodeNumber, &this->splitV[nodeNumber], &this->splitP[nodeNumber], this->csplit, NULL, this->nodes[rightChild] ,
                                         this->data, this->nInstances, this->nVariables,   this->variables);
    }else if( rightChild <= 0){
      this->nodes[nodeNumber] = new Node(nodeNumber, &this->splitV[nodeNumber], &this->splitP[nodeNumber], this->csplit,this->nodes[leftChild], NULL,
                                         this->data, this->nInstances, this->nVariables,  this->variables);
    }else{
      this->nodes[nodeNumber] = new Node(nodeNumber, &this->splitV[nodeNumber], &this->splitP[nodeNumber],this->csplit, this->nodes[leftChild], this->nodes[rightChild],
                                         this->data, this->nInstances, this->nVariables, this->variables);
    }
  }else{
    this->nodes[nodeNumber] = NULL;
  }
  
} // end initNode


int Tree::predictClass(int minBucket, int minSplit, bool pruneIfInvalid, int nodeNumber){
  // calculate predictions
  // if pruneIfInvalid == TRUE non-valid nodes a pruned
  // otherwise -1 is returned for non-valid nodes
  if(nodeNumber == 0){
    for(int i = 0; i < *this->nInstances; i++){
      this->classification[i] = 0;
    }
  }else{
    this->reverseClassification(nodeNumber, nodeNumber);
  }
  
  int returnValue = this->nodes[nodeNumber]->partition( this->classification, this->weights, this->variables, &this->nNodes, minBucket, minSplit);
  if(returnValue == -1){
    return -1;
  }else if(returnValue <= 0 || pruneIfInvalid == false){
    return returnValue;
  }else{
    this->deleteChildNodes(returnValue); // call recursion delete node and everything below it
    return predictClass(minBucket, minSplit, true, 0);
  }
} // end  predictClass


bool Tree::reverseClassification(int startNode, int nodeNumber){
  // all observations which are in the subtree below "startNode" are classified to be in node "startNode"
  // this saves computation time when the split-rule in startNode is changed. After the call of "reverseClassification" only the instances in
  // node "startNode" are newly evaluated. The rest of the tree stays the same.
  for(int i = 0; i < *this->nInstances; i++){
    if(this->classification[i] == (nodeNumber*2+1) || this->classification[i] == (nodeNumber*2+2) ){
      this->classification[i] = startNode;
    }
  }
  if(nodeNumber*2+1 < *this->maxNode){
    if(splitV[nodeNumber] >= 0){
      reverseClassification(startNode, nodeNumber*2+1);
    }
  }
  if(nodeNumber*2+2 < *this->maxNode){
    if(splitV[nodeNumber] >= 0){
      reverseClassification(startNode, nodeNumber*2+2);
    }
  }
  return true;
} // end reverseClassification


bool Tree::deleteChildNodes(int nodeNumber){
  // used by predictClass() and mutateNode() to delete a child node
  if(this->splitV[nodeNumber] >= 0  && nodeNumber > 0){
    if(this->nodes[nodeNumber]->leftChild != NULL ){
      this->deleteChildNodes(nodeNumber*2+1);
    }
    if(this->nodes[nodeNumber]->rightChild != NULL ){
      this->deleteChildNodes(nodeNumber*2+2);
    }
    if(nodeNumber%2 == 0){
      this->nodes[(int) ((nodeNumber-1)/2)]->rightChild = NULL;
    }else{
      this->nodes[(int) ((nodeNumber-1)/2)]->leftChild = NULL;
    }
    this->splitV[nodeNumber] = -999999;
    this->splitP[nodeNumber] = -999999;
    this->nNodes--;
    delete this->nodes[nodeNumber];
    this->nodes[nodeNumber] = NULL;
    return true;
  }else{
    // cout << "warning: node could not be deleted " << endl;
    return false;
  }
} // end deleteChildNodes


void Tree::randomizeCategories(int nodeNumber){
  // assigns random categories of a categorical variable
  bool left = false;
  bool right = false;
  for(int i = 0; i < this->variables[ *this->nodes[nodeNumber]->splitV ]->nCats ; i++){
    if(i == this->variables[*this->nodes[nodeNumber]->splitV ]->nCats-1 && left == false){
      this->csplit[i][nodeNumber] = 1;
    }else if(i == this->variables[*this->nodes[nodeNumber]->splitV]->nCats-1 && right == false){
      this->csplit[i][nodeNumber] = 3;
    }else if(getUnifRandNumber(2) == 1){
      this->csplit[i][nodeNumber] = 1;
      left = true;
    }else{
      this->csplit[i][nodeNumber] = 3;
      right = true;
    }
  }
  
} // end randomizeCategories


bool Tree::calculateTotalCosts(int method, double alpha, double beta, double clv, double d, double f, double lambda, int sumWeights){
  // calculates tree quality
  IntegerVector nObsClass1PerNode;
  IntegerVector nObsPerNode;
  calculateAllScores(0, nObsClass1PerNode, nObsPerNode);
  
  std::vector<double> scores_std;
  scores_std.reserve(sumWeights);
  std::vector<int> classes_std;
  classes_std.reserve(sumWeights);
  
  for(int i = 0; i < nObsPerNode.size(); i++){
    double scorePerNode_i = nObsClass1PerNode[i] / (double) nObsPerNode[i];
    NumericVector vect_score_in_node(nObsPerNode[i], scorePerNode_i);
    scores_std.insert(scores_std.end(), vect_score_in_node.begin(), vect_score_in_node.end());
    
    IntegerVector vect_class1_in_node(nObsClass1PerNode[i], 1);
    IntegerVector vect_class0_in_node(nObsPerNode[i] - nObsClass1PerNode[i], 0);
    classes_std.insert(classes_std.end(), vect_class1_in_node.begin(), vect_class1_in_node.end());
    classes_std.insert(classes_std.end(), vect_class0_in_node.begin(), vect_class0_in_node.end());
  }
  
  NumericVector scores = wrap(scores_std);
  IntegerVector classes = wrap(classes_std);
  
  this->performance = -empChurnCpp(scores, classes, alpha, beta, clv, d, f) + lambda*(this->nNodes+1.0);
  
  return true;
} // end calculateTotalCosts


void Tree::calculateAllScores(int nodeNumber, IntegerVector& nObsClass1PerNode, IntegerVector& nObsPerNode){
  // calculates all scores and reports all corresponding classes
  if(this->nodes[nodeNumber]->leftChild != NULL){
    this->calculateAllScores(nodeNumber*2+1, nObsClass1PerNode, nObsPerNode);
  }
  if(this->nodes[nodeNumber]->rightChild != NULL){
    this->calculateAllScores(nodeNumber*2+2, nObsClass1PerNode, nObsPerNode);
  }
  
  if( this->splitV[nodeNumber] >= 0 && this->nodes[nodeNumber]->leftChild == NULL){
    this->nodes[nodeNumber]->calculateChildNodeScore(true, this->weights, nObsClass1PerNode, nObsPerNode);
  }
  if( this->splitV[nodeNumber] >= 0 && this->nodes[nodeNumber]->rightChild == NULL){
    this->nodes[nodeNumber]->calculateChildNodeScore(false, this->weights, nObsClass1PerNode, nObsPerNode);
  }
  
} // end calculateAllScores

