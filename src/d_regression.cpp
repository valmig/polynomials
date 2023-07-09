#include <pol.h>
#include <string>
#include <val_utils.h>
#include <d_array.h>
#include <LA.h>
#include <fstream>
#include <pol_arithmetic.h>


void program_info()
{
	std::cout<<"\nCommand: d_regression Input-File/\"pairs of numbers\" [deg = 1]";
	std::cout<<"\n=======";
	std::cout<<"\nComputes the regression double-polynomial of degree deg";
	std::cout<<"\nwrt the points given by the pairs of numbers!\n";
}



val::d_array<double> GetPoints(const std::string& s)
{
	val::d_array<double> Points;
	
	int n=s.length();
	if (!n) return Points;
	std::string sn;
	
	val::d_array<double>::set_plus_cap(10);
	
	
	for (int i=0;i<n;++i) {
		if (s[i]==' ' || s[i]== '\n' || s[i]==';') {
			Points.push_back(val::FromString<double>(sn));
			sn="";
		}
		else sn+=s[i];
	}
	Points.push_back(val::FromString<double>(sn));
		
	return Points; 
}

val::pol<double> regression(const val::d_array<double>& Points,int degree=1)
{
	val::pol<double> f;
	int n = Points.length(), N=n/2,i,j;
	
	if (n<2) return f;
	degree = val::Min(degree,N);
    
    val::d_array<val::vector<double> > g(degree+1);
    val::vector<double> y(N);
    double det;

    for (i=0;i<=degree;++i) g[i]=val::vector<double>(N);
    for (i=0;i<N;++i) g[0](i)=1.0;
    for (i=0;i<N;++i) {
         g[1](i) = Points[2*i];
         y(i) = Points[2*i+1];
    }
	for (i=2;i<=degree;++i) {
		for (j=0;j<N;++j) g[i](j) = val::power(g[1](j),i);
	}

	val::matrix<double> A(degree+1,degree+2),X;
	for (i=0;i<=degree;++i) {
		A(i,i) = g[i]*g[i];
		for (j=i+1;j<=degree;++j) A(i,j) = A(j,i) = g[i]*g[j];
		A(i,degree+1) = g[i]*y;
	}
	val::les_double(A,X,det);

    for (i=0;i<=degree;++i) f.insert(X(0,i),i);
    return f;
}


int main(int argnr,char* argv[])
{
    //test_interpolation();

	std::string firstarg;
	if (argnr == 2) firstarg = std::string(argv[1]);
	
	if (argnr == 1 || (argnr == 2 && (firstarg == "-h" || firstarg == "--help") )) {
        program_info();
        return 0;
    }
    
	std::string s;
	std::fstream file(argv[1], std::ios::in);
	
	if (file) { 
		std::string line;
		while (std::getline(file,line)) {
			s += line + " ";
		}
	}
	else s = std::string(argv[1]);
	

    val::d_array<double> Points = GetPoints(s);
    int degree=1;
    if (argnr>2) degree = val::FromString<int>(std::string(argv[2]));
    if (degree<1) degree=1;
    
    val::pol<double> f = regression(Points,degree);
    std::cout<<"\nRegression polynomial:\n"<<f<<"\n = "<<val::PolToString(f)<<"\n\n";
        
    return 0;
}
