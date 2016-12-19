#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include "payoff.h"
#include "option_value.h"
#include "binomial.h"
#include "trinomial.h"

int main(int argc, char* argv[])
{
	// parameters
	double T = 8.0 / 12.0;
	double K = 48;
	double S = 50;
	double r = 0.02;
	double q = 0.01;
	double vol = 0.3;
	int N = 8;

	// print precision
	int p = 9;

	CallPayoff barrierCall(K);
	BinomialTree b1(barrierCall, T, S, r, q, vol);



	

	std::ofstream ofs;
	/*
	ofs.open("res.txt");
	ofs << std::fixed
	<< std::setprecision(p) << b1.barrierOption_binomial_tree(1000, 45).price;
	
	for(int i = 10; i<=1000 ;i++){
		ofs << b1.barrierOption_binomial_tree(i, 45).price
			<< std::endl;
	}
	*/


	//b1.barrierOption_exact(K, T, S, r, q, vol, 45);
	
	ofs.open("res_tri.txt");
	std::cout << std::fixed
		<< std::setprecision(p);
	/*
	TrinomialTree t1(barrierCall, T, S, r, q, vol);
	for (int i = 10; i <= 1000; i++) {
		ofs << t1.barrierOption_trinomial_tree(i,45).price
			<< std::endl;
	}
	*/



	//find optimal N for binomial tree
	std::vector<double> K_value;
	int size = 0;
	for (int i = 10; i <= 1000; i++)
	{
		K_value.push_back(int(std::log(45/ S) / (vol*sqrt(T / i))) + 1);
		size++;
		if (size > 1)
		{
			if (K_value[size - 1] == K_value[size - 2] - 1)
			{
				std::cout << i - 1 << ","  << b1.barrierOption_binomial_tree(i-1, 45).price << std::endl;
			}
		}
	}
	
	



}