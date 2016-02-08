#include <constant_parameters.h>

Constant_parameters::Constant_parameters(double constant)
                                        : d_constant(constant),
                                          d_constant_square(constant*constant)
{}

double Constant_parameters::constant() const
{
    return d_constant;
}

double Constant_parameters::constant_square() const
{
    return d_constant_square;
}

Parameters* Constant_parameters::clone() const
{
    return new Constant_parameters(*this);
}

double Constant_parameters::integral(double time1, double time2) const
{
    return d_constant*(time2-time1);
}

double Constant_parameters::square_integral(double time1, double time2) const
{
    return d_constant_square*(time2-time1);
}

double Constant_parameters::time_product_diff(double time1,
                                              double time2) const
{
    return d_constant*(time2-time1);
}
