#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <valarray>
#include <fstream>
#include <string>

using namespace std;

class ode1{
public:
    double get_t0(){return t0;}
    double get_tn(){return tn;}
    double get_h(){return h;}
    double get_n(){return n;}
    valarray<double> get_x(){return x;}
    string get_filename(){return filename;}
    // derivative is redefined in derived classes
    virtual valarray<double> derivative(const double& t, const valarray<double>& x){return {};}
protected:
    // domain and partitions
    double t0; // start time
    double tn; // end time
    double h; // step size
    int n; // number of steps
    // initial condition
    valarray<double> x;
    // file name for output
    string filename;
    // default constructor / destructor is protected and cannot be used in main()
    ode1(){}
    ~ode1(){}
};

class lorenz : public ode1{
public:
    // lorenz(string filename = "lorenz") : filename(filename){} // does not work
    // lorenz(string filename = "lorenz"){filename = filename;} does not work
    lorenz(string fn = "lorenz"){
        filename = fn;
        setupParameter();
        setupDomainAndPartition();
        setupInitialCondition();
    }
    ~lorenz(){}
    // valarray<double> derivative(valarray<double>& x){} does not work because (w + k2) is not an lvalue
    valarray<double> derivative(const double& t, const valarray<double>& x){
        return {sigma * (x[1] - x[0]), R * x[0] - x[1] - x[0] * x[2], x[0] * x[1] - b * x[2]};
    }
private:
    // parameters
    double sigma;
    double R;
    double b;
    // setup
    void setupParameter();
    void setupDomainAndPartition();
    void setupInitialCondition();
};

class cooling : public ode1{
public:
    cooling(string fn = "cooling"){
        filename = fn;
        setupParameter();
        setupDomainAndPartition();
        setupInitialCondition();
    }
    ~cooling(){}
    valarray<double> derivative(const double& t, const valarray<double>& T){
        return {k * (Ts - T)};
    }
private:
    // parameters
    double k; // heat transfer coefficient
    double Ts; // surrounding temperature
    // setup
    void setupParameter();
    void setupDomainAndPartition();
    void setupInitialCondition();
};

class RungeKutta{
public:
    RungeKutta(int order = 4) : order(order){}
    ~RungeKutta(){}
    // void solve(const lorenz& lorenz){} does not work
    void solve(ode1& sys){
        file.open(sys.get_filename() + "-rk" + to_string(order) + ".dat", std::ios::in|std::ios::out|std::ios::trunc);
        ti = sys.get_t0();
        w = sys.get_x();
        output(file, ti, w);
        if(order == 1){
            // Forward Euler
            for(int i = 0; i < sys.get_n(); ++i){
                k1 = sys.get_h() * sys.derivative(ti, w);
                w = w + k1;
                ti = ti + sys.get_h();
                output(file, ti, w);
            }
        }else if(order == 2){
            // Modified Euler Method
            for(int i = 0; i < sys.get_n(); ++i){
                k1 = sys.get_h() * sys.derivative(ti, w);
                k2 = sys.get_h() * sys.derivative(ti + sys.get_h(), w + k1);
                w = w + (k1 + k2) / 2;
                ti = ti + sys.get_h();
                output(file, ti, w);
            }
        }else if(order == 3){
            // Heun's method
            for(int i = 0; i < sys.get_n(); ++i){
                k1 = sys.get_h() * sys.derivative(ti, w);
                k2 = sys.get_h() * sys.derivative(ti + sys.get_h() / 3, w + k1 / 3);
                k3 = sys.get_h() * sys.derivative(ti + sys.get_h() * 2 / 3, w + k2 * 2 / 3);
                w = w + (k1 + 3 * k3) / 4;
                ti = ti + sys.get_h();
                output(file, ti, w);
            }
        }else{
            // RK4
            for(int i = 0; i < sys.get_n(); ++i){
                k1 = sys.get_h() * sys.derivative(ti, w);
                k2 = sys.get_h() * sys.derivative(ti + sys.get_h() / 2, w + k1 / 2);
                k3 = sys.get_h() * sys.derivative(ti + sys.get_h() / 2, w + k2 / 2);
                k4 = sys.get_h() * sys.derivative(ti + sys.get_h(), w + k3);
                w = w + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                ti = ti + sys.get_h();
                output(file, ti, w);
            }
        }
        file.close();
    }
private:
    int order;
    double ti; // current time
    valarray<double> w;
    valarray<double> k1;
    valarray<double> k2;
    valarray<double> k3;
    valarray<double> k4;
    fstream file;
    void output(fstream& file, const double& ti, const valarray<double>& w){
        file << ti;
        for(int i = 0; i < w.size(); ++i){
            file << " " << w[i];
        }
        file << endl;
    }
};

#endif // ODE_SOLVER_H