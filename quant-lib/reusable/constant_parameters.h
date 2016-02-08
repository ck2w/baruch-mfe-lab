#ifndef CONSTANT_PARAMETERS_H
#define CONSTANT_PARAMETERS_H
#include <parameters.h>

class Constant_parameters : public Parameters
{
    public:

        Constant_parameters(double constant);

        double constant() const;
        double constant_square() const;

        virtual Parameters* clone() const;
        virtual double integral(double time1, double time2) const;
        virtual double square_integral(double time1, double time2) const;
        virtual double time_product_diff(double time1, double time2) const;

    private:

        double d_constant;
        double d_constant_square;
};

#endif /* CONSTANT_PARAMETERS_H */
