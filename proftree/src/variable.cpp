#include "variable.h"

variable::variable(int varNumber, int dependendVar, int nInst, double** data, int nCats){
  // init Variables
  if(nCats < 0){
    this->isCat = true;
    this->nCats= -nCats;
  }else{
    this->isCat = false;
    this->nCats = nCats;
  }
  this->sortedValues = new double[this->nCats ];
  for(int i = 0; i < this->nCats; i++){
    this->sortedValues[i]= -999999;
  }
  if(varNumber != dependendVar && this->isCat == false){
    this->sortedValues[0] = data[0][varNumber];
    int i = 1;
    int k = 1;
    while( i < nInst && k < this->nCats){
      int n = 0;
      bool flag = false;
      while(n < k && flag ==  false){
        if(data[i][varNumber] ==  this->sortedValues[n]){
          flag = true;
        }
        n++;
      }
      if(flag == false){
        this->sortedValues[k] = data[i][varNumber] ;
        k++;
      }
      i++;
    }
    this->sortValues();
  }else if(varNumber != dependendVar){
    for(int i = 0; i < this->nCats; i++)
      this->sortedValues[i] = i+1;
  }
} // end variable


variable::~variable(){
  // deconstructor destroys the instance and frees up memory
  delete [] sortedValues;
  sortedValues = NULL;
} // end ~variable


void variable::sortValues(){
  // sort the values of a variable
  double* temp = new double[this->nCats];
  for( int i = 0; i < this->nCats; i++){
    int count = 0;
    for( int j = 0; j < this->nCats; j++){
      if(this->sortedValues[i] > this->sortedValues[j])
        count ++;
    }
    temp[count] = this->sortedValues[i];
  }
  for( int i = 0; i < this->nCats; i++){
    this->sortedValues[i] = temp[i];
  }
  delete [] temp;
  temp= NULL;
} // end sortValues

