#ifndef CUMULATIVE_NORMAL_H
#define CUMULATIVE_NORMAL_H

class Cumulative_normal
{
    public:

        static double approximate_evaluate(double x);
        static double exact_evaluate(double x);

    private:

        static const double a0 = 0.2316419;
        static const double a1 = 0.319381530;
        static const double a2 = -0.356563782;
        static const double a3 = 1.781477937;
        static const double a4 = -1.821255978;
        static const double a5 = 1.330274429;
};

#endif /* CUMULATIVE_NORMAL_H */
