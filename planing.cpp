#include <math.h>
#include "planning.h"
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


//路径规划入口函数
double Path::CalcP(Pinput& pinput){
    double dx = pinput.x1 - pinput.x0;
	double dy = pinput.y1[0] - pinput.y0[0];

    b.clear();
	b.push_back(pinput.y0[0]);
	b.push_back(pinput.y0[1]);
	b.push_back(pinput.y0[2] / 2);
	double A2[3][3] = { { pow(dx, 3), pow(dx, 4), pow(dx, 5) },
	{ 3 * pow(dx, 2), 4 * pow(dx, 3), 5 * pow(dx, 4) },
	{ 6 * dx, 12 * pow(dx, 2), 20 * pow(dx, 3) } };
	double B2[3] = { pinput.y1[0] - b[0] - b[1] * dx - b[2] * pow(dx, 2), pinput.y1[1] - b[1] - 2 * b[2] * dx, pinput.y1[2] - 2 * b[2] };
	double C2[3] = { 0 };
	double A2_inversion[3][3] = { 0 };
	Matrix_inversion(A2, A2_inversion);
	Matrix_multiple(A2_inversion, B2, C2);
	b.push_back(C2[0]);
	b.push_back(C2[1]);
	b.push_back(C2[2]);

    double s = S_to_X(pinput);
    
	return s;
}



//建立s-x对应表，插值时通过s找到x
double Path::S_to_X(Pinput& pinput){
	double stepx = 0.1;
    double dx = pinput.x1 - pinput.x0;
	int m = ceil(dx / stepx) + 1;
	double *x, *y;
	x = (double*)malloc(sizeof(double)* m);
	y = (double*)malloc(sizeof(double)* m);
	double temp = 0;
	x[0] = pinput.x0;
	y[0] = b[0];
    table_x.clear();
    table_s.clear();
    table_x.push_back(x[0]);
    table_s.push_back(0);
	for (int i = 1; i < m; i++)
	{
		x[i] = x[0] + stepx*i;
		if (x[i]>pinput.x1)
		{
			x[i] = pinput.x1;
		}
		y[i] = b[0] + b[1] * x[i] + b[2] * pow(x[i], 2) + b[3] * pow(x[i], 3) + b[4] * pow(x[i], 4) + b[5] * pow(x[i], 5);
		temp = sqrt(pow(x[i] - x[i - 1], 2) + pow(y[i] - y[i - 1], 2));
        table_x.push_back(x[i]);
        table_s.push_back(table_s[i-1] + temp);
	}

    return table_s[m-1];
}

/*
函数功能：根据s值查找x；
参数：s
返回：s处的x
*/
double Path::Xr(double s){
    int n = 0;
    for(int i=0; i<table_s.size()-1; i++){
        if(s>=table_s[i] && s<table_s[i+1]){
            n = i;
            break;
        }
    }
    double xr = table_x[n] + (s - table_s[n])/(table_s[n+1] - table_s[n])*(table_x[n+1] - table_x[n]); 
    return xr;
}


/*
函数功能：根据x值求y；
参数：x
返回：x处的y
*/
double Path::Yr(double x){
    double yr = b[0] + b[1] * x + b[2] * pow(x, 2) + b[3] * pow(x, 3) + b[4] * pow(x, 4) + b[5] * pow(x, 5);
    return yr;
}

/*
函数功能：根据x值求theta；
参数：x
返回：x处的theta
*/
double Path::Thetar(double x){
	double yr_dx = b[1] + 2 * b[2] * x + 3 * b[3] * pow(x, 2) + 4 * b[4] * pow(x, 3) + 5 * b[5] * pow(x, 4);
	double thetar = atan(yr_dx);
    return thetar;
}

/*
函数功能：根据x值求曲率k；
参数：x
返回：x处的k
*/
double Path::Kr(double x){
    double yr_dx = b[1] + 2 * b[2] * x + 3 * b[3] * pow(x, 2) + 4 * b[4] * pow(x, 3) + 5 * b[5] * pow(x, 4);
    double yr_dx2 = 2 * b[2] + 6 * b[3] * x + 12 * b[4] * pow(x, 2) + 20 * b[5] * pow(x, 3);
	double kr = yr_dx2 / pow(1 + pow(yr_dx, 2), 1.5);
    return kr;
}


/*
函数功能：求矩阵的逆
参数1：原矩阵
参数2：逆矩阵
*/
void Path::Matrix_inversion(double mat[3][3], double mat_inversion[3][3])
{
	int i, j;
	double determinant = 0;
	//求行列式
	for (i = 0; i < 3; i++)
		determinant = determinant + (mat[0][i] * (mat[1][(i + 1) % 3] * mat[2][(i + 2) % 3] - mat[1][(i + 2) % 3] * mat[2][(i + 1) % 3]));
	//求逆
	for (i = 0; i < 3; i++){
		for (j = 0; j < 3; j++)
			mat_inversion[i][j] = ((mat[(j + 1) % 3][(i + 1) % 3] * mat[(j + 2) % 3][(i + 2) % 3]) - (mat[(j + 1) % 3][(i + 2) % 3] * mat[(j + 2) % 3][(i + 1) % 3])) / determinant;
	}
}


/*
函数功能：3x3矩阵乘3x1矩阵
参数：a*b=c
*/
void Path::Matrix_multiple(double a[3][3], double b[3], double c[3])
{
	int m = 3, n = 3;
	for (int i = 0; i <= m - 1; i++)
	{
		c[i] = 0;
		for (int j = 0; j <= n - 1; j++)
		{
			c[i] = a[i][j] * b[j] + c[i];
		}
	}
}



//速度规划入口函数
bool Speed::CalcS(Param& param){
	if (!InitParam(param)){
		return false;
	}

	param.vmin = -param.vmax;
	param.amin = -param.amax;
	param.jmin = -param.jmax;
	double T1 = sqrt(abs(param.v1 - param.v0) / param.jmax);
	double T2 = param.amax / param.jmax;
	double Tjs = min(T1, T2);
	double Dq = param.q1 - param.q0;
	double Tj1;
	double Tj2;
	double Ta;
	double Td;
	double Tv;
	double Tj;
	double Delta;
	if (T1 <= T2){
		if (Dq < Tjs*(param.v0 + param.v1)){
			cout << "位移过小，不存在满足始末速度的轨迹！";
			return false;
		}
	}
	else{
		if (Dq < 0.5*(param.v0 + param.v1)*(Tjs + abs(param.v1 - param.v0) / param.amax)){
			cout << "位移过小，不存在满足始末速度的轨迹！";
		    return false;
        }
	}
	//输入参数正确误（即轨迹存在），分类讨论
	if ((param.vmax - param.v0)*param.jmax < pow(param.amax, 2)){         //(3.19)满足,amax不能达到
		Tj1 = sqrt((param.vmax - param.v0) / param.jmax);
		Ta = 2 * Tj1;
		param.alima = param.jmax*Tj1;
	}
	else{       //(3.19)不满足，amax能达到
		Tj1 = param.amax / param.jmax;
		Ta = Tj1 + (param.vmax - param.v0) / param.amax;
		param.alima = param.amax;
	}

	if ((param.vmax - param.v1)*param.jmax < pow(param.amax, 2)){       //(3.20)满足,amin不能达到
		Tj2 = sqrt((param.vmax - param.v1) / param.jmax);
		Td = 2 * Tj2;
		param.alimd = -param.jmax*Tj2;
	}
	else{                                                       //(3.20)不满足,amin能达到
		Tj2 = param.amax / param.jmax;
		Td = Tj2 + (param.vmax - param.v1) / param.amax;
		param.alimd = -param.amax;
	}

	//计算（3.25）
	Tv = (param.q1 - param.q0) / param.vmax - Ta / 2 * (1 + param.v0 / param.vmax) -
		Td / 2 * (1 + param.v1 / param.vmax);
	if (Tv > 0){                                              //case1,最大速度能达到
		param.vlim = param.vmax;
		param.Tj1 = Tj1;
		param.Ta = Ta;
		param.Tj2 = Tj2;
		param.Td = Td;
		param.T = Ta + Tv + Td;
		param.Tv = Tv;
		return true;
	}
	else{                                                      // case2,最大速度不能达到
		Tv = 0;
		param.Tv = Tv;
		//计算（3.26a）,(3.27),（3.26b）,（3.26c）
		Tj = param.amax / param.jmax;
		Tj1 = Tj;
		Tj2 = Tj;
		Delta = pow(param.amax, 4) / pow(param.jmax, 2) + 2 * (pow(param.v0, 2) + pow(param.v1, 2)) + param.amax*
			(4 * (param.q1 - param.q0) - 2 * param.amax / param.jmax*(param.v0 + param.v1));
		Ta = (pow(param.amax, 2) / param.jmax - 2 * param.v0 + sqrt(Delta)) / (2 * param.amax);
		Td = (pow(param.amax, 2) / param.jmax - 2 * param.v1 + sqrt(Delta)) / (2 * param.amax);
		if (Ta > 2 * Tj && Td > 2 * Tj){                                 //加速段和减速段都能达到最大加速度
			param.Tj1 = Tj1;
			param.Tj2 = Tj2;
			param.Ta = Ta;
			param.Td = Td;
			param.T = Ta + Tv + Td;
			param.alima = param.amax;
			param.alimd = -param.amax;
			param.vlim = param.v0 + (Ta - Tj1)*param.alima;
			return true;
		}
		else{                                                        //至少有一段不能达到最大加速度
			double gamma = 0.99;
			double amax = param.amax;
			//逐渐减小最大加速度约束
			while (Ta < 2 * Tj || Td < 2 * Tj){
				if (Ta > 0 && Td > 0){
					amax = gamma*amax;
					//循环计算（3.26a）,(3.27),（3.26b）,（3.26c）
					Tj = amax / param.jmax;
					Tj1 = Tj;
					Tj2 = Tj;
					Delta = pow(amax, 4) / pow(param.jmax, 2) + 2 * (pow(param.v0, 2) + pow(param.v1, 2)) + amax*
						(4 * (param.q1 - param.q0) - 2 * amax / param.jmax*(param.v0 + param.v1));
					Ta = (pow(amax, 2) / param.jmax - 2 * param.v0 + sqrt(Delta)) / (2 * amax);
					Td = (pow(amax, 2) / param.jmax - 2 * param.v1 + sqrt(Delta)) / (2 * amax);
				}
				else{                                                  //出现Ta或Td小于0
					if (Ta <= 0){
						Ta = 0;
						Tj1 = 0;
						//计算（3.28a）
						Td = 2 * (param.q1 - param.q0) / (param.v0 + param.v1);
						//计算（3.28b）
						double num = param.jmax*(param.q1 - param.q0) - sqrt(param.jmax*(param.jmax*pow((param.q1 - param.q0), 2) +
							pow((param.v1 + param.v0), 2)*(param.v1 - param.v0)));
						double den = param.jmax*(param.v1 + param.v0);
						Tj2 = num / den;
					}
					else if (Td <= 0){
						Td = 0;
						Tj2 = 0;
						//计算（3.29a）
						Ta = 2 * (param.q1 - param.q0) / (param.v0 + param.v1);
						//计算（3.29b）
						double num = param.jmax*(param.q1 - param.q0) - sqrt(param.jmax*(param.jmax*pow((param.q1 - param.q0), 2) -
							pow((param.v1 + param.v0), 2)*(param.v1 - param.v0)));
						double den = param.jmax*(param.v1 + param.v0);
						Tj1 = num / den;
					}
					param.Tj1 = Tj1;
					param.Tj2 = Tj2;
					param.Ta = Ta;
					param.Td = Td;
					param.T = Ta + Tv + Td;
					param.alima = param.jmax*Tj1;
					param.alimd = -param.jmax*Tj2;
					param.vlim = param.v0 + (Ta - Tj1)*param.alima;
					return true;
				}
			}
			param.Tj1 = Tj1;
			param.Tj2 = Tj2;
			param.Ta = Ta;
			param.Td = Td;
			param.T = Ta + Tv + Td;
			param.alima = param.jmax*Tj1;
			param.alimd = -param.jmax*Tj2;
			param.vlim = param.v0 + (Ta - Tj1)* param.alima;
			return true;
		}
	}
}

bool Speed::InitParam(Param& param){
	//参数检查，绝对值小于1.0E-8的数值认为是0.
	if (abs(param.vmax) < 1.0E-8){
		cout << "参数错误。最大速度不能为0";
		return false;
	}
	if (abs(param.amax) < 1.0E-8){
		cout << "参数错误。最大加速度不能为0";
		return false;
	}
	if (abs(param.jmax) < 1.0E-8){
		cout << "参数错误。最大加加速度不能为0";
		return false;
	}
	if (abs(param.q1 - param.q0) < 1.0E-8){
		cout << "参数错误。输入位移不能为0";
		return false;
	}
	return true;
}

//返回时间t对应的s
double Speed::Dis(double t, Param& param){
	double q = 0;
	if (t >= 0 && t < param.Tj1) {
		q = param.q0 + param.v0 * t + param.jmax * pow(t, 3) / 6;
	}
	else if (t >= param.Tj1 && t < (param.Ta - param.Tj1)) {
		q = param.q0 + param.v0 * t + param.alima / 6 * (3 * pow(t, 2) - 3 * param.Tj1 * t + pow(param.Tj1, 2));
	}
	else if (t >= param.Ta - param.Tj1 && t < param.Ta) {
		q = param.q0 + (param.vlim + param.v0) * param.Ta / 2 - param.vlim * (param.Ta - t) - param.jmin * pow((param.Ta - t), 3) / 6;
	}
	else if (t >= param.Ta && t < param.Ta + param.Tv) {
		q = param.q0 + (param.vlim + param.v0) * param.Ta / 2 + param.vlim * (t - param.Ta);
	}
	else if (t >= param.T - param.Td && t < param.T - param.Td + param.Tj2) {
		q = param.q1 - (param.vlim + param.v1) * param.Td / 2 + param.vlim * (t - param.T + param.Td) - param.jmax * pow((t - param.T + param.Td), 3) / 6;
	}
	else if (t >= param.T - param.Td + param.Tj2 && t < param.T - param.Tj2) {
		q = param.q1 - (param.vlim + param.v1) * param.Td / 2 + param.vlim * (t - param.T + param.Td) + param.alimd / 6 * (3 * pow((t - param.T + param.Td), 2) - 3 * param.Tj2 * (t - param.T + param.Td) + pow(param.Tj2, 2));
	}
	else {
		q = param.q1 - param.v1 * (param.T - t) - param.jmax * pow((param.T - t), 3) / 6;
	}
	return q;
}

//返回时间t对应的v
double Speed::Velocity(double t, Param& param){
	double v0 = param.v0;
	double v1 = param.v1;
	double vlim = param.vlim;
	double alima = param.alima;
	double alimd = param.alimd;
	double jmax = param.jmax;
	double Tj1 = param.Tj1;
	double Ta = param.Ta;
	double Tv = param.Tv;
	double Tj2 = param.Tj2;
	double Td = param.Td;
	double T = param.T;
	double jmin = param.jmin;
	double dq = 0;
	if (t >= 0 && t < Tj1){
		dq = v0 + jmax*pow(t, 2) / 2;
	}
	else if (t >= Tj1 && t < Ta - Tj1){
		dq = v0 + alima*(t - Tj1 / 2);
	}
	else if (t >= Ta - Tj1 && t < Ta){
		dq = vlim + jmin*pow((Ta - t), 2) / 2;
	}
	else if (t >= Ta && t < Ta + Tv){
		dq = vlim;
	}
	else if (t >= T - Td && t < T - Td + Tj2){
		dq = vlim - jmax*pow((t - T + Td), 2) / 2;
	}
	else if (t >= T - Td + Tj2 && t < T - Tj2){
		dq = vlim + alimd*(t - T + Td - Tj2 / 2);
	}
	else{
		dq = v1 + jmax*pow((T - t), 2) / 2;
	}
	return dq;
}

//返回时间t对应的a
double Speed::Acceleration(double t, Param& param){
	double alima = param.alima;
	double alimd = param.alimd;
	double jmax = param.jmax;
	double Tj1 = param.Tj1;
	double Ta = param.Ta;
	double Tv = param.Tv;
	double Tj2 = param.Tj2;
	double Td = param.Td;
	double T = param.T;
	double jmin = param.jmin;
	double ddq = 0;
	if (t >= 0 && t < Tj1){
		ddq = jmax*t;
	}
	else if (t >= Tj1 && t < Ta - Tj1){
		ddq = alima;
	}
	else if (t >= Ta - Tj1 && t < Ta){
		ddq = -jmin*(Ta - t);
	}
	else if (t >= Ta && t < Ta + Tv){
		ddq = 0;
	}
	else if (t >= T - Td && t < T - Td + Tj2){
		ddq = -jmax*(t - T + Td);
	}
	else if (t >= T - Td + Tj2 && t < T - Tj2){
		ddq = alimd;
	}
	else{
		ddq = -jmax*(T - t);
	}
	return ddq;
}

//返回时间t对应的j
double Speed::Jerk(double t, Param& param){
	double jmax = param.jmax;
	double Tj1 = param.Tj1;
	double Ta = param.Ta;
	double Tv = param.Tv;
	double Tj2 = param.Tj2;
	double Td = param.Td;
	double T = param.T;
	double jmin = param.jmin;
	double jt = 0;
	if (t >= 0 && t < Tj1){
		jt = jmax;
	}
	else if (t >= Tj1 && t < Ta - Tj1){
		jt = 0;
	}
	else if (t >= Ta - Tj1 && t < Ta){
		jt = jmin;
	}
	else if (t >= Ta && t < Ta + Tv){
		jt = 0;
	}
	else if (t >= T - Td && t < T - Td + Tj2){
		jt = jmin;
	}
	else if (t >= T - Td + Tj2 && t < T - Tj2){
		jt = 0;
	}
	else{
		jt = jmax;
	}
	return jt;
}




//输入：
// xstart   起点x坐标
// xend     终点x坐标
// ystart   起点y位置，切线斜率(可默认为0)，曲率相关(默认为0)
// yend     终点y位置，切线斜率(可默认为0)，曲率相关(默认为0)
// v0  v1   始末速度
// amax     最大加速度
// vmax     最大速度
// jmax     最大加加速度

// 输出：
// time     时间
// dis      路程s
// vel      速度
// acc      加速度
// jerk     加加速度
// xr       x
// yr       y
// thetar   theta
// kr       曲率

int main(){
    double xstart = 0;
	double xend = 100;
	double ystart[3] = { 0, 0, 0 };
	double yend[3] = { 0, 0, 0 };

	time.clear();
	dis.clear();
	vel.clear();
	acc.clear();
	jerk.clear();
	xr.clear();
	yr.clear();
	thetar.clear();
	kr.clear();

    Path P1;
    Pinput pinput = {xstart, xend, { 0, 0, 0 }, { 0, 0, 0 }};
    double s = P1.CalcP(pinput);
	cout<<"s = "<<s<<endl;
	Param param = {0};
	param.q0 = 0;
	param.q1 = s;
	param.v0 = 0;
	param.v1 = 15;
	param.vmax = 15;
	param.amax = 4;
	param.jmax = 1;
    Speed S1;
	if (!S1.CalcS(param)){
		cout << "planning failed";
	}
	cout<<"T = "<< param.T<<endl;

	double t = 0;
	int i = 0;
	while(t < param.T){
		 time.push_back(t);
	     dis.push_back(S1.Dis(t,param));
		 vel.push_back(S1.Velocity(t, param));
		 acc.push_back(S1.Acceleration(t, param));
		 jerk.push_back(S1.Jerk(t, param));
         xr.push_back(P1.Xr(dis[i]));
         yr.push_back(P1.Yr(dis[i]));
         thetar.push_back(P1.Thetar(dis[i]));
         kr.push_back(P1.Kr(dis[i]));
         if(vel[i]<0.1){
             t=t+1;
         }else if(0.1<=vel[i] && vel[i]<2){
             t=t+0.5;
         }else if(2<=vel[i] && vel[i]<5){
             t=t+0.2;
         }else if(5<=vel[i] && vel[i]<10){
             t=t+0.15;
         }else{
             t=t+0.1;
         }
         i++;
	 }
	time.push_back(param.T);
	dis.push_back(S1.Dis(param.T,param));
	vel.push_back(S1.Velocity(param.T, param));
	acc.push_back(S1.Acceleration(param.T, param));
	jerk.push_back(S1.Jerk(param.T, param));
	xr.push_back(P1.Xr(s));
	yr.push_back(P1.Yr(s));
	thetar.push_back(P1.Thetar(s));
	kr.push_back(P1.Kr(s));
	for(int j=0; j<dis.size(); j++){	
		 cout<< j <<"   ";
		 cout<<"t = "<< time[j] <<"   ";
		 cout<<"s = "<< dis[j] <<"   ";
		 cout<<"v = "<< vel[j] <<"  ";
		 cout<<"a = "<< acc[j] <<"  ";
		 cout<<"j = "<< jerk[j] <<"  ";
		 cout<<"xr = "<< xr[j] <<"  ";
		 cout<<"yr = "<< yr[j] <<"  ";
		 cout<<"thetar = "<< thetar[j] <<"  ";
		 cout<<"kr = "<< kr[j] << endl;
	}

	getchar();
	return 0;
}

