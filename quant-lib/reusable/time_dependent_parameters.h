#ifndef TIME_DEPENDENT_PARAMETERS_H
#define TIME_DEPENDENT_PARAMETERS_H
#include <parameters.h>

class Time_dependent_parameters : public Parameters
{
    public:

        Time_dependent_parameters( double (*func)(double) );

        virtual Parameters* clone() const;
        virtual double integral(double time1, double time2) const;
        virtual double square_integral(double time1, double time2) const;
        virtual double time_product_diff(double time1, double time2) const;

    private:

        double (*d_func)(double);
};

#endif /* TIME_DEPENDENT_PARAMETERS_H */
