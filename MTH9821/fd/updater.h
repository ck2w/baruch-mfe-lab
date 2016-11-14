#ifndef UPDATER_H
#define UPDATER_H 
#include <vector>
#include <Eigen/Dense>

class Updater
{
    public:

        Updater() : d_alpha(0) {}

        double getAlpha() const
        {
            return d_alpha;
        }
        
        void setAlpha(double c)
        {
            d_alpha = c;
        }

        virtual void config(double c, int N)=0;

        virtual void update( const std::vector<double> & uOld,
                             std::vector<double>* uNew ) const=0;

    private:

        double d_alpha;
};

class ForwardEulerUpdater : public Updater
{
    public:

        virtual void config(double c, int N);

        virtual void update( const std::vector<double> & uOld,
                             std::vector<double>* uNew ) const;
};

class BackwardEulerUpdater : public Updater
{
    public:

        virtual void config(double c, int N);

        virtual void update( const std::vector<double> & uOld,
                             std::vector<double>* uNew ) const;

    private:

        Eigen::MatrixXd d_L;
        Eigen::MatrixXd d_U;
};

class BackwardEulerSorUpdater : public Updater
{
    public:

        BackwardEulerSorUpdater( double w ) : d_omega(w), d_tol(1e-6) {}

        virtual void config(double c, int N);

        virtual void update( const std::vector<double> & uOld,
                             std::vector<double>* uNew ) const;

    private:

        double d_omega;
        double d_tol;

        Eigen::MatrixXd d_A;
};

class CrankNicolsonUpdater : public Updater
{
    public:

        virtual void config(double c, int N);

        virtual void update( const std::vector<double> & uOld,
                             std::vector<double>* uNew ) const;

    private:

        Eigen::MatrixXd d_L;
        Eigen::MatrixXd d_U;
        Eigen::MatrixXd d_B;
};

class CrankNicolsonSorUpdater : public Updater
{
    public:
        
        CrankNicolsonSorUpdater( double w ) : d_omega(w), d_tol(1e-6) {}

        virtual void config(double c, int N);

        virtual void update( const std::vector<double> & uOld,
                             std::vector<double>* uNew ) const;

    private:

        double d_omega;
        double d_tol;

        Eigen::MatrixXd d_A;
        Eigen::MatrixXd d_B;
};

#endif /* UPDATER_H */
