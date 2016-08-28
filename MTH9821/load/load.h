#ifndef LOAD_H
#define LOAD_H 

#include <Eigen/Dense>
#include <string>

int countRows(const char * filename);

Eigen::MatrixXd readCSV(const char* file, int rows, int cols);

Eigen::VectorXd readCSV(const char* file, int size);

#endif /* LOAD_H */
