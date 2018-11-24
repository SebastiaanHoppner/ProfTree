#include "node.h"

Node::Node(int pos, int* splitV, double* splitP, int** csplit, Node* leftChild, Node* rightChild, double** data, int* nInst, int* nVar, variable** variables){
  this->pos = pos;
  this->splitV = splitV;
  this->splitP = splitP;
  this->nInst = nInst;
  this->nVar = nVar;
  this->data = data;
  this->leftChild = leftChild;
  this->rightChild = rightChild;
  this->variables = variables;
  this->localClassification = new int[*this->nInst];
  for(int i = 0; i < *this->nInst; i++)
    this->localClassification[i] = 0;
  this->sumLocalWeights = 0;
  this->sumLeftLocalWeights = 0;
  this->sumRightLocalWeights = 0;
  this->predictionInternalNode = 0;
  this->predictionLeftTerminal = 0;
  this->predictionRightTerminal = 0;
  this->nClassesDependendVar = &variables[*this->nVar-1]->nCats;
  this->csplit = csplit;
} // end Node


int Node::partition( int* classification, int* weights, variable** variables, int* nNodes, int minBucket, int minSplit ){
  // assigns instances to belong to the right or the left child node
  for(int i = 0; i < *this->nInst; i++)
    this->localClassification[i] = classification[i];
  this->sumLeftLocalWeights = 0;
  this->sumRightLocalWeights = 0;
  if(this->variables[*this->splitV]->isCat == true){ // categorical split-variable
    bool flag = false;
    for(int i = 0; i < *this->nInst; i++){
      if(classification[i] == this->pos){
        flag = false;
        for(int k = 0; k < variables[ *this->splitV]->nCats && flag == false ; k++){
          if( variables[*this->splitV ]->sortedValues[k] == this->data[i][*this->splitV] ){
            if( this->csplit[k][this->pos] == 1 ){
              classification[i] = (this->pos)*2+1;
              this->localClassification[i] = classification[i];
              this->sumLeftLocalWeights++;
            }else{
              classification[i] = (this->pos)*2+2;
              this->localClassification[i] = classification[i];
              this->sumRightLocalWeights++;
            }
            flag = true;
          }
        }
      }
    }
  }else if( variables[*this->splitV]->isCat == false){  // numeric split-variable
    for(int i = 0; i < *this->nInst; i++){
      if(classification[i] == this->pos){
        if( (double)this->data[i][*this->splitV ] < (double)*this->splitP){
          classification[i] = (this->pos)*2+1;
          this->sumLeftLocalWeights += weights[i];
          this->localClassification[i] = classification[i];
        }else{
          classification[i] = (this->pos)*2+2;
          this->sumRightLocalWeights += weights[i];
          this->localClassification[i] = classification[i];
        }
      }
    }
    
  }
  this->sumLocalWeights = this->sumLeftLocalWeights + this->sumRightLocalWeights ;
  if( this->sumLocalWeights < minSplit  && this->pos > 0){  // checks if there are enough instances in the nodes; otherwise the node is pruned
    return (int) this->pos;
  }
  
  // recursive call of partition until there are no further internal nodes
  int temp1 = -1, temp2 = -1;
  if( this->leftChild != NULL){
    temp1 = this->leftChild->partition(classification, weights, variables, nNodes, minBucket, minSplit);
  }
  if( this->rightChild != NULL){
    temp2 = this->rightChild->partition(classification, weights, variables, nNodes, minBucket, minSplit);
  }
  
  if( temp1 == -2 || temp2 == -2){
    return -2;
  }else if( temp1 == 0 || temp2 == 0){
    return 0;
  }else if(temp1 != -1){
    return temp1;
  }else if(temp2 != -1){
    return temp2;
  }
  
  if(this->sumLeftLocalWeights < minBucket){   // if there are to few instances in the left terminal-node the internal node is pruned
    return this->pos;
  }else if(this->sumRightLocalWeights < minBucket){ // if there are to few instances in the right terminal-node the internal node is pruned
    return this->pos;
  }
  return -1;
} // end partition


void Node::calculateChildNodeScore(bool leftNode, int* weights, IntegerVector& nObsClass1PerNode, IntegerVector& nObsPerNode){
  // calculates the total number of observations AND the number of class-1 observations in a node
  int sumWeights = 0; // sum of weights of all observations in the child node
  int sumOnes = 0; // sum of weights of class-1 observations in the child node
  
  if(leftNode == true){
    for(int i = 0; i < *this->nInst; i++){
      if(localClassification[i] == (this->pos)*2+1){
        sumWeights += weights[i];
        // std::cout << (int) data[i][*this->nVar-1] << "\n";
        if ((int) data[i][*this->nVar-1] == 2)
          sumOnes += weights[i];
      }
    }
  }else{
    for(int i = 0; i < *this->nInst; i++){
      if(localClassification[i] == (this->pos)*2+2 ){
        sumWeights += weights[i];
        // std::cout << (int) data[i][*this->nVar-1] << "\n";
        if ((int) data[i][*this->nVar-1] == 2)
          sumOnes += weights[i];
      }
    }
  }
  
  nObsClass1PerNode.push_back(sumOnes);
  nObsPerNode.push_back(sumWeights);
} // end calculateChildNodeScore


Node::~Node(){
  // deconstructor destroys the instance and frees up memory
  delete [] localClassification;
  localClassification = NULL;
  leftChild = NULL;
  rightChild = NULL;
  splitV = NULL;
  splitP = NULL;
  csplit = NULL;
  leftChild = NULL;
  rightChild = NULL;
  nInst = NULL;
  nVar = NULL;
  data = NULL;
} // ~end Node

