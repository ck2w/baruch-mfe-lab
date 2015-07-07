#ifndef PARAMETERS_H
#define PARAMETERS_H 

class Parameters
{
    public:

        virtual Parameters* clone() const=0;
        virtual double integral(double time1, double time2) const=0;
        virtual double square_integral(double time1, double time2) const=0;
        virtual double time_product_diff(double time1, double time2) const=0;

        virtual ~Parameters() {}
};

class Parameters_bridge
{
    public:

        Parameters_bridge( const Parameters_bridge & parameters_bridge );
        Parameters_bridge( const Parameters & parameters );

        Parameters_bridge& operator=( const Parameters_bridge & parameters_bridge );

        double integral(double time1, double time2) const;
        double square_integral(double time1, double time2) const;
        double time_product_diff(double time1, double time2) const;

        double mean(double time1, double time2) const;
        double rms(double time1, double time2) const;

        ~Parameters_bridge();

    private:

        Parameters* d_parameters_ptr;
};

#endif /* PARAMETERS_H */
