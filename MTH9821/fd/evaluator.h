#ifndef EVALUATOR_H
#define EVALUATOR_H 

class Evaluator
{
    public:

        Evaluator() {}

        virtual double operator()(double s) const
        {
            return 0.0;
        }

        virtual double operator()(double x, double t) const
        {
            return 0.0;
        }
};

class FunctionEvaluator : public Evaluator
{
    public:

        FunctionEvaluator(double (*f)(double)) : d_f(f) {}

        virtual double operator()(double s) const
        {
            return (*d_f)(s);
        }
        
    private:

        double (*d_f)(double);
};

#endif /* EVALUATOR_H */
