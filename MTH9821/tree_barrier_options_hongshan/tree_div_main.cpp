#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include "payoff.h"
#include "option_value.h"
#include "binomial.h"

int main(int argc, char* argv[])
{
	// parameters
	double T = 8.0 / 12.0;
	double K = 55.55;
	double S = 50;
	double r = 0.02;
	double q = 0;
	double vol = 0.3;

	double prev;//record previous value of before twice the N
	double cur;
	int N; //initial numeer of steps
	
	int p = 9;// print precision
	double tol = 1e-4;


	
	//(1)

	//Three Month European Put :
	PutPayoff EuroPut(55.55);
	BinomialTree p_3month_euro(EuroPut, 3.0/12.0, S, r, q, vol);
	N = 12;
	prev = p_3month_euro.evaluate_prop_div(6, 0.01, 2.0 / 12.0).price;
	cur = p_3month_euro.evaluate_prop_div(N, 0.01, 2.0 / 12.0).price;
	while (std::abs(cur - prev) >= tol) {
		prev = cur;
		N *= 2;
		cur = p_3month_euro.evaluate_prop_div(N, 0.01, 2.0 / 12.0).price;

	}
	
	auto res = p_3month_euro.evaluate_prop_div(N, 0.01, 2.0 / 12.0);
	std::cout << "Three Month European Put:"
		<< std::fixed
		<< std::setprecision(p)
		<< std::endl
		<< res.price
		<< "," << N
		<< "," << res.delta
		<< "," << res.gamma
		<< "," << res.theta
		<< std::endl;
	

	/*
	//Three Month American put :
	N = 12;
	prev = p_3month_euro.evaluate_prop_div(6, 0.01, 2.0 / 12.0,true).price;
	cur = p_3month_euro.evaluate_prop_div(N, 0.01, 2.0 / 12.0,true).price;
	while (std::abs(cur - prev) >= tol) {
		prev = cur;
		N *= 2;
		res = p_3month_euro.evaluate_prop_div(N, 0.01, 2.0 / 12.0, true);
		cur = res.price;
		
	}

	std::cout << "Three Month American Put:"
		<< std::fixed
		<< std::setprecision(p)
		<< std::endl
		<< res.price
		<< "," << N
		<< "," << res.delta
		<< "," << res.gamma
		<< "," << res.theta
		<< std::endl;

	*/

	
	//Seven Euro put :
	std::cout << "seven Month American Put:";
	BinomialTree p_7month_euro(EuroPut, 7.0 / 12.0, S, r, q, vol);
	N = 28672;
	//prev = p_7month_euro.evaluate_prop_div(7, 0.01, 2.0 / 12.0).price;
	res = p_7month_euro.evaluate_fix_div(N, 0.5, 2.0 / 12.0, true);
	//res = p_7month_euro.evaluate_fix_prop_div_mix(N,true);
	

	/*
	while (std::abs(cur - prev) >= tol) {
		
		prev = cur;
		N *= 2;
		res = p_3month_euro.evaluate_prop_div(N, 0.01, 2.0 / 12.0);
		cur = res.price;
		std::cout << std::abs(cur - prev)  << std::endl;

	}
	*/

	/*
	//7 month Euro put fixed div :

	BinomialTree p_3month_euro_fix(EuroPut, 3.0 / 12.0, S, r, q, vol);
	N = 192;
	res = p_3month_euro_fix.evaluate_fix_div(N,0.5, 2.0 / 12.0);
	*/

	std::ofstream ofs;
	ofs.open("p3_res.txt", std::ofstream::app);
	ofs << std::fixed
		<< std::setprecision(p);
	
	ofs << "7 Month american fix dividend:"
	<< std::fixed
	<< std::setprecision(p)
	<< std::endl
	<< res.price
	<< "," << N
	<< "," << res.delta
	<< "," << res.gamma
	<< "," << res.theta
	<< std::endl;

	


}