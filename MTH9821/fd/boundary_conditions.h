#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H 
#include <evaluator.h>

class VanillaCallTerminalCondition : public Evaluator
{
    public:

        VanillaCallTerminalCondition() {}
        VanillaCallTerminalCondition( double expiry,
                                      double spot,
                                      double strike,
                                      double rate,
                                      double div,
                                      double vol )
        {
            d_strike = strike;
            double c = (rate-div)/(vol*vol);
            d_a = c - 0.5;
            d_b = (c+0.5)*(c+0.5) + 2*div/(vol*vol);
        }

        virtual double operator()(double s) const
        {
            return d_strike*std::exp(d_a*s)*std::max(std::exp(s)-1,0.0);
        }
        
        virtual void setBoundaries(double ti, double tf, double xl, double xr) 
        {}
        
    private:

        double d_strike;
        double d_a;
        double d_b;
};

typedef VanillaCallTerminalCondition AmericanCallTerminalCondition;

class VanillaCallRightBoundaryCondition : public Evaluator
{
    public:

        VanillaCallRightBoundaryCondition() {}
        VanillaCallRightBoundaryCondition( double expiry,
                                           double spot,
                                           double strike,
                                           double rate,
                                           double div,
                                           double vol )
                                         : d_expiry(expiry),
                                           d_spot(spot),
                                           d_strike(strike),
                                           d_rate(rate),
                                           d_div(div),
                                           d_vol(vol)
        {
            double c = (rate-div)/(vol*vol);
            d_a = c - 0.5;
            d_b = (c+0.5)*(c+0.5) + 2*div/(vol*vol);
        }

        virtual double operator()(double s) const
        {
            double c1 = std::exp(-2*d_rate*s/(d_vol*d_vol));
            double c2 = std::exp(d_xr - 2*d_div*s/(d_vol*d_vol));
            return d_strike*std::exp(d_a*d_xr+d_b*s)*(c2-c1); 
        }
        
        virtual void setBoundaries(double ti, double tf, double xl, double xr) 
        {
            d_xr = xr; 
        }

    private:

        double d_expiry;
        double d_spot;
        double d_strike;
        double d_rate;
        double d_div;
        double d_vol;
        double d_a;
        double d_b;
        double d_xr;
};

class VanillaCallLeftBoundaryCondition : public Evaluator
{
    public:

        VanillaCallLeftBoundaryCondition() {}
        VanillaCallLeftBoundaryCondition( double expiry,
                                          double spot,
                                          double strike,
                                          double rate,
                                          double div,
                                          double vol )
        {}

        virtual double operator()(double s) const
        {
            return 0; 
        }
        
        virtual void setBoundaries(double ti, double tf, double xl, double xr) 
        {} 
};

typedef VanillaCallLeftBoundaryCondition AmericanCallLeftBoundaryCondition;

class AmericanCallRightBoundaryCondition : public Evaluator
{
    public:

        AmericanCallRightBoundaryCondition() {}
        AmericanCallRightBoundaryCondition( double expiry,
                                            double spot,
                                            double strike,
                                            double rate,
                                            double div,
                                            double vol )
                                          : d_expiry(expiry),
                                            d_spot(spot),
                                            d_strike(strike),
                                            d_rate(rate),
                                            d_div(div),
                                            d_vol(vol)
        {
            double c = (rate-div)/(vol*vol);
            d_a = c - 0.5;
            d_b = (c+0.5)*(c+0.5) + 2*div/(vol*vol);
        }

        virtual double operator()(double s) const
        {
            return d_strike*d_c*std::exp(d_a*d_xr+d_b*s);
        }
        
        virtual void setBoundaries(double ti, double tf, double xl, double xr) 
        {
            d_xr = xr; 
            d_c = std::exp(d_xr)-1;
        }

    private:

        double d_expiry;
        double d_spot;
        double d_strike;
        double d_rate;
        double d_div;
        double d_vol;
        double d_a;
        double d_b;
        double d_c;
        double d_xr;
};

class VanillaPutTerminalCondition : public Evaluator
{
    public:

        VanillaPutTerminalCondition() {}
        VanillaPutTerminalCondition( double expiry,
                                     double spot,
                                     double strike,
                                     double rate,
                                     double div,
                                     double vol )
        {
            d_strike = strike;
            double c = (rate-div)/(vol*vol);
            d_a = c - 0.5;
            d_b = (c+0.5)*(c+0.5) + 2*div/(vol*vol);
        }

        virtual double operator()(double s) const
        {
            return d_strike*std::exp(d_a*s)*std::max(1-std::exp(s),0.0);
        }
        
        virtual void setBoundaries(double ti, double tf, double xl, double xr) 
        {} 
        
    private:

        double d_strike;
        double d_a;
        double d_b;
};

typedef VanillaPutTerminalCondition AmericanPutTerminalCondition;

class VanillaPutLeftBoundaryCondition : public Evaluator
{
    public:

        VanillaPutLeftBoundaryCondition() {}
        VanillaPutLeftBoundaryCondition( double expiry,
                                         double spot,
                                         double strike,
                                         double rate,
                                         double div,
                                         double vol )
                                       : d_expiry(expiry),
                                         d_spot(spot),
                                         d_strike(strike),
                                         d_rate(rate),
                                         d_div(div),
                                         d_vol(vol)
        {
            double c = (rate-div)/(vol*vol);
            d_a = c - 0.5;
            d_b = (c+0.5)*(c+0.5) + 2*div/(vol*vol);
        }

        virtual double operator()(double s) const
        {
            double c1 = std::exp(-2*d_rate*s/(d_vol*d_vol));
            double c2 = std::exp(d_xl - 2*d_div*s/(d_vol*d_vol));
            return d_strike*std::exp(d_a*d_xl+d_b*s)*(c1-c2); 
        }
        
        virtual void setBoundaries(double ti, double tf, double xl, double xr) 
        {
            d_xl = xl;
        } 
        
    private:

        double d_expiry;
        double d_spot;
        double d_strike;
        double d_rate;
        double d_div;
        double d_vol;
        double d_a;
        double d_b;
        double d_xl;
};

class VanillaPutRightBoundaryCondition : public Evaluator
{
    public:

        VanillaPutRightBoundaryCondition() {}
        VanillaPutRightBoundaryCondition( double expiry,
                                          double spot,
                                          double strike,
                                          double rate,
                                          double div,
                                          double vol )
        {}

        virtual double operator()(double s) const
        {
            return 0; 
        }
        
        virtual void setBoundaries(double ti, double tf, double xl, double xr) 
        {}
};

typedef VanillaPutRightBoundaryCondition AmericanPutRightBoundaryCondition;

class AmericanPutLeftBoundaryCondition : public Evaluator
{
    public:

        AmericanPutLeftBoundaryCondition() {}
        AmericanPutLeftBoundaryCondition( double expiry,
                                          double spot,
                                          double strike,
                                          double rate,
                                          double div,
                                          double vol )
                                        : d_expiry(expiry),
                                          d_spot(spot),
                                          d_strike(strike),
                                          d_rate(rate),
                                          d_div(div),
                                          d_vol(vol)
        {
            double c = (rate-div)/(vol*vol);
            d_a = c - 0.5;
            d_b = (c+0.5)*(c+0.5) + 2*div/(vol*vol);
        }

        virtual double operator()(double s) const
        {
            return d_strike*d_c*std::exp(d_a*d_xl+d_b*s);
        }
        
        virtual void setBoundaries(double ti, double tf, double xl, double xr) 
        {
            d_xl = xl;
            d_c = 1-std::exp(d_xl);
        }
        
    private:

        double d_expiry;
        double d_spot;
        double d_strike;
        double d_rate;
        double d_div;
        double d_vol;
        double d_a;
        double d_b;
        double d_c;
        double d_xl;
};

class AmericanPutEarlyExercisePremium : public Evaluator
{
    public:

        AmericanPutEarlyExercisePremium() {}
        AmericanPutEarlyExercisePremium( double expiry,
                                         double spot,
                                         double strike,
                                         double rate,
                                         double div,
                                         double vol )
                                       : d_expiry(expiry),
                                         d_spot(spot),
                                         d_strike(strike),
                                         d_rate(rate),
                                         d_div(div),
                                         d_vol(vol)
        {
            double c = (rate-div)/(vol*vol);
            d_a = c - 0.5;
            d_b = (c+0.5)*(c+0.5) + 2*div/(vol*vol);
        }

        virtual double operator()(double x, double t) const
        {
            return d_strike*std::exp(d_a*x+d_b*t)*std::max(1-std::exp(x), 0.0);
        }

        virtual void setBoundaries(double ti, double tf, double xl, double xr) {}
        
    private:

        double d_expiry;
        double d_spot;
        double d_strike;
        double d_rate;
        double d_div;
        double d_vol;
        double d_a;
        double d_b;
};

class VanillaCallWithDividendTerminalCondition : public Evaluator
{
    public:

        VanillaCallWithDividendTerminalCondition() {}
        VanillaCallWithDividendTerminalCondition( double expiry,
                                                  double spot,
                                                  double strike,
                                                  double rate,
                                                  double div,
                                                  double vol,
                                                  double qDiv )
        {
            d_strike = strike;
            double c = (rate-div)/(vol*vol);
            d_a = c - 0.5;
            d_b = (c+0.5)*(c+0.5) + 2*div/(vol*vol);
            d_shift = -std::log(1-qDiv);
        }

        virtual double operator()(double s) const
        {
            return d_strike*std::exp(d_a*(s+d_shift))*std::max(std::exp(s)-1,0.0);
        }
        
        virtual void setBoundaries(double ti, double tf, double xl, double xr) 
        {}
        
    private:

        double d_strike;
        double d_a;
        double d_b;
        double d_shift;
};

class VanillaCallWithDividendRightBoundaryCondition : public Evaluator
{
    public:

        VanillaCallWithDividendRightBoundaryCondition() {}
        VanillaCallWithDividendRightBoundaryCondition( double expiry,
                                                       double spot,
                                                       double strike,
                                                       double rate,
                                                       double div,
                                                       double vol,
                                                       double qDiv )
                                                     : d_expiry(expiry),
                                                       d_spot(spot),
                                                       d_strike(strike),
                                                       d_rate(rate),
                                                       d_div(div),
                                                       d_vol(vol)
        {
            d_strike = strike;
            double c = (rate-div)/(vol*vol);
            d_a = c - 0.5;
            d_b = (c+0.5)*(c+0.5) + 2*div/(vol*vol);
            d_shift = -std::log(1-qDiv);
        }

        virtual double operator()(double s) const
        {
            double c1 = std::exp(-2*d_rate*s/(d_vol*d_vol));
            double c2 = std::exp(d_xr - 2*d_div*s/(d_vol*d_vol));
            return d_strike*std::exp(d_a*(d_xr+d_shift)+d_b*s)*(c2-c1); 
        }
        
        virtual void setBoundaries(double ti, double tf, double xl, double xr) 
        {
            d_xr = xr; 
        }
        
    private:

        double d_expiry;
        double d_spot;
        double d_strike;
        double d_rate;
        double d_div;
        double d_vol;
        double d_a;
        double d_b;
        double d_xr;
        double d_shift;
};

class VanillaCallWithDividendLeftBoundaryCondition : public Evaluator
{
    public:

        VanillaCallWithDividendLeftBoundaryCondition() {}
        VanillaCallWithDividendLeftBoundaryCondition( double expiry,
                                                      double spot,
                                                      double strike,
                                                      double rate,
                                                      double div,
                                                      double vol,
                                                      double qDiv )
        {}

        virtual double operator()(double s) const
        {
            return 0; 
        }
        
        virtual void setBoundaries(double ti, double tf, double xl, double xr) 
        {} 
};

#endif /* BOUNDARY_CONDITIONS_H */
