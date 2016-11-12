#include <ran.h>
#include <gtest/gtest.h>
#include <iostream>
#include <cmath>
#include <vector>

class RanTest : public ::testing::Test
{
    protected:
        
        virtual void SetUp() {}
        virtual void TearDown() {}
};

TEST_F(RanTest, PrintSomeRandomNumbers)
{
    int idum=1;
    int count=0;
    int ns = 10000;
    for (int i=0; i<ns; i++) {
        double num = bmnormal(&idum, &count);
    }
}

/*
TEST_F(RanTest, CheckRandomNumberGenerator)
{
    int idum=1;
    int ns = 10000000;
    int ni = 100;
    
    double num=0.0;
    double interval = 1.0/ni;
    std::vector<int> count(ni,0);

    for (int i=0; i<ns; i++) {
        num = ran(&idum);
        int index = num/interval;
        count[index]++;
    }

    double maxDeviation = 0.0;
    for (int i=0; i<ni; i++) {
        double density = count[i]/double(ns);
        double deviation = std::fabs(density - interval);
        if (deviation > maxDeviation) {
            maxDeviation = deviation;
        }
    }

    double expectedDeviation = 1.0/std::sqrt(ns);
    EXPECT_LT(maxDeviation, expectedDeviation);
}

TEST_F(RanTest, CheckInverseTransformStandardNormal)
{
    int idum=1;
    int count=0;
    int ns = 10000000;
    int ni = 100;
    double interval = 0.05;
    
    double num=0.0;
    std::vector<int> countPos(ni,0);
    std::vector<int> countNeg(ni,0);
    int countZero = 0;

    for (int i=0; i<ns; i++) {
        num = normal(&idum, &count);
        int index = num/interval;
        if (index>0) {
            if (index>=ni) {
                index = ni;
            }
            countPos[index-1]++;
        }
        else if (index<0) {
            index = -index;
            if (index >= ni) {
                index = ni;
            }
            countNeg[index-1]++;
        }
        else {
            countZero++;
        }
    }

    EXPECT_NEAR(double(countZero)/(countPos[0]+countNeg[0]), 1.0, 0.01);
}

TEST_F(RanTest, CheckAcceptanceRejectionStandardNormal)
{
    int idum=1;
    int count=0;
    int ns = 10000000;
    int ni = 100;
    double interval = 0.05;
    
    double num=0.0;
    std::vector<int> countPos(ni,0);
    std::vector<int> countNeg(ni,0);
    int countZero = 0;

    for (int i=0; i<ns; i++) {
        num = arnormal(&idum, &count);
        int index = num/interval;
        if (index>0) {
            if (index>=ni) {
                index = ni;
            }
            countPos[index-1]++;
        }
        else if (index<0) {
            index = -index;
            if (index >= ni) {
                index = ni;
            }
            countNeg[index-1]++;
        }
        else {
            countZero++;
        }
    }

    EXPECT_NEAR(double(countZero)/(countPos[0]+countNeg[0]), 1.0, 0.01);
}

TEST_F(RanTest, CheckBoxMullerStandardNormal)
{
    int idum=1;
    int count=0;
    int ns = 200;
    int ni = 100;
    double interval = 0.05;
    
    double num=0.0;
    std::vector<int> countPos(ni,0);
    std::vector<int> countNeg(ni,0);
    int countZero = 0;

    for (int i=0; i<ns; i++) {
        num = bmnormal(&idum, &count);
        int index = num/interval;
        if (index>0) {
            if (index>=ni) {
                index = ni;
            }
            countPos[index-1]++;
        }
        else if (index<0) {
            index = -index;
            if (index >= ni) {
                index = ni;
            }
            countNeg[index-1]++;
        }
        else {
            countZero++;
        }
    }

    EXPECT_NEAR(double(countZero)/(countPos[0]+countNeg[0]), 1.0, 0.05);
}
*/
