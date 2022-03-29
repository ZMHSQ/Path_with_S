#include <vector>


struct Pinput{
    double x0;
    double x1;
    double y0[3];
    double y1[3];
};


class Path
{
public:
    double CalcP(Pinput& pinput);
    double S_to_X(Pinput& pinput);
    double Xr(double s);
    double Yr(double x);
    double Thetar(double x);
    double Kr(double x);
    void Matrix_inversion(double mat[3][3], double mat_inversion[3][3]);
    void Matrix_multiple(double a[3][3], double b[3], double c[3]);

private:
    std::vector<double> b;
    // std::vector<double> xr;
    // std::vector<double> yr;
    // std::vector<double> thetar;
    // std::vector<double> kr;
    std::vector<double> table_s;
    std::vector<double> table_x;
};


struct Param{
	double Tj1;
	double Tj2;
	double Ta;
	double Tv;
	double Td;
	double T;
	double q0;
	double q1;
	double v0;
	double v1;
	double vmax;
	double vmin;
	double amax;
	double amin;
	double jmax;
	double jmin;
	double vlim;
	double alima;
	double alimd;

	// Param(){}
	// Param(double Tj1_, double Tj2_, double Ta_, double Tv_, double Td_, double T_,
	//         double q0_, double q1_, double v0_, double v1_, double vlim_, double alima_, 
	//         double alimd_, double jmax_, double jmin_) : 
	//         Tj1(Tj1_), Tj2(Tj2_), Ta(Ta_), Tv(Tv_), Td(Td_), T(T_), q0(q0_), q1(q1_), 
	//         v0(v0_), v1(v1_), vlim(vlim_), alima(alima_), alimd(alimd_), jmax(jmax_), jmin(jmin_) {}


};

class Speed
{
public:
	bool CalcS(Param& param);

	bool InitParam(Param& param);

	double Dis(double t, Param& param);
	double Velocity(double t, Param& param);
	double Acceleration(double t, Param& param);
	double Jerk(double t, Param& param);

private:
	// std::vector<double> dis;
	// std::vector<double> vel;
	// std::vector<double> acc;
	// std::vector<double> jerk;

};
	std::vector<double> time;	
	std::vector<double> dis;
	std::vector<double> vel;
	std::vector<double> acc;
	std::vector<double> jerk;
    std::vector<double> xr;
    std::vector<double> yr;
    std::vector<double> thetar;
    std::vector<double> kr;