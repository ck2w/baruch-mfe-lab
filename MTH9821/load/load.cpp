#include <load.h>
#include <Eigen/Dense>
#include <fstream>

int countRows(const char * filename)
{
    int numRows = 0;
    std::string unused;
    std::ifstream in(filename);
    while (std::getline(in, unused)) {
        if (!unused.empty()) {
            numRows++;
        }
    }

    return numRows;
}

Eigen::MatrixXd readCSV(const char* file, int rows, int cols)
{
    std::ifstream in(file);
    std::string line;

    int row = 0;
    int col = 0;

    Eigen::MatrixXd res = Eigen::MatrixXd(rows, cols);
    if (in.is_open()) {
        while (std::getline(in, line)) {
            char *ptr = (char *) line.c_str();
            int len = line.length();
            col = 0;
            char *start = ptr;
            for (int i=0; i<len; i++) {
                if (ptr[i] == ',') {
                    res(row, col++) = atof(start);
                    start = ptr + i + 1;
                }
            }
            res(row, col) = atof(start);
            row++;
        }

        in.close();
    }

    return res;
}

Eigen::VectorXd readCSV(const char* file, int size)
{
    std::ifstream in(file);
    std::string line;

    Eigen::VectorXd res = Eigen::VectorXd(size);
    if (in.is_open()) {
        std::getline(in, line);
        char *ptr = (char *) line.c_str();
        int len = line.length();
        int col = 0;
        char *start = ptr;
        for (int i=0; i<len; i++) {
            if (ptr[i] == ',') {
                res(col++) = atof(start);
                start = ptr + i + 1;
            }
        }
        res(col) = atof(start);
        in.close();
    }

    return res;
}

