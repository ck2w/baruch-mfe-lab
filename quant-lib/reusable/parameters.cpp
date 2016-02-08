#include <parameters.h>

Parameters_bridge::Parameters_bridge( const Parameters_bridge & parameters_bridge )
{
    d_parameters_ptr = parameters_bridge.d_parameters_ptr->clone();
}

Parameters_bridge::Parameters_bridge( const Parameters & parameters )
{
    d_parameters_ptr = parameters.clone();
}

Parameters_bridge& Parameters_bridge::operator=( const Parameters_bridge & parameters_bridge )
{
    if (this != &parameters_bridge)
    {
        delete d_parameters_ptr;
        d_parameters_ptr = parameters_bridge.d_parameters_ptr->clone();
    }

    return *this;
}

Parameters_bridge::~Parameters_bridge()
{
    delete d_parameters_ptr;
}

double Parameters_bridge::integral(double time1, double time2) const
{
    return d_parameters_ptr->integral(time1, time2);
}

double Parameters_bridge::square_integral(double time1, double time2) const
{
    return d_parameters_ptr->square_integral(time1, time2);
}

double Parameters_bridge::time_product_diff(double time1, double time2) const
{
    return d_parameters_ptr->time_product_diff(time1, time2);
}

double Parameters_bridge::mean(double time1, double time2) const
{
    return integral(time1, time2)/(time2-time1);
}

double Parameters_bridge::rms(double time1, double time2) const
{
    return square_integral(time1, time2)/(time2-time1);
}
