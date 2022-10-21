#ifndef GODUNOV_SOLVER_GODUNOV_SOLVER_HPP
#define GODUNOV_SOLVER_GODUNOV_SOLVER_HPP

class solver {
public:
    solver() = default;

    ~solver() = default;

    void solve_system(double x_0, double x_n, double x_c, double time, int N,
                      double r, double u, double p);

    void solve_system(double x_0, double x_n, double x_c, double time, int N,
                      double r_l, double u_l, double p_l,
                      double r_r, double u_r, double p_r);

private:
    void output_parameters_to_file();
};

#endif //GODUNOV_SOLVER_GODUNOV_SOLVER_HPP
