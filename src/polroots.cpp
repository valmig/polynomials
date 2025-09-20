#include <fstream>
#include <analysis.h>
#include <vector.h>
#include <complex.h>
#include <polfactor.h>
#include <valfunction.h>
#include <d_array.h>
#include <pol_arithmetic.h>


void program_info()
{
	std::cout<<"\nCommand: polroots [-r] \"polynomial\"/Input-File";
	std::cout<<"\n=======";
	std::cout<<"\nTries to compute an approximation of all real and ";
	std::cout<<"\ncomplex roots of a polynomial with real coefficients!";
	std::cout<<"\nWith the -r option rational roots are computed!\n";
}


int isinfix(const std::string &s)
{
    int n=s.size(),i;

    for (i=0;i<n;i++)
        if (s[i]=='x') return 1;
    return 0;

}


int main(int argnr,char* argv[])
{
	using namespace val;
	
	std::string firstarg;
	int computerational = 0;
	
	if (argnr >= 2) firstarg = std::string(argv[1]);
	if (argnr == 1 || (argnr == 2 && (firstarg == "-h" || firstarg == "--help") )) {
        program_info();
        return 0;
	}
	else if (argnr == 3) {
		if (firstarg == "-r") {
			computerational = 1;
			firstarg = std::string(argv[2]);
		}
		else {
			std::cout << "\n Invalid argument: " << firstarg;
			return 1;
		}
	}		 
	else if (argnr > 3) {
		std::cout << "\n Too many arguments. Quit program!!\n";
		return 1;
	}
 
 	std::string s,line;
	val::pol<val::rational> f_r;
	
	if (firstarg.length() <= 20) {
		std::fstream file(firstarg, std::ios::in);
		if (file) {
			while (std::getline(file,line)) {
				s += line + " ";
			}
		}
		else s = firstarg;
	}
	else  s= firstarg;
	
	//std::cout << "\n s = "<<s;
	
	if (isinfix(s)) {
		val::valfunction F(s);
		f_r = F.getpolynomial();
	}
	else f_r = val::FromString<val::pol<val::rational>>(s);


	val::d_array<val::rational> Rationalzeros;

	if (computerational) Rationalzeros= rational_roots(f_r);     // declared in polfactor.h; defined in polfactor.cpp
	Rationalzeros.sort();
 
 	val::pol<double> f = val::ToDoublePolynom(f_r);

	val::vector<double> Realzeros;
	val::vector<val::complex> Compzeros;
	std::string os;
	int i,NR,NC,N;
 
 	//val::realRoots(f,Realzeros);
	//std::cout<<"\nReal zeros computed!"<<std::endl;

	N=val::roots(f,Realzeros,Compzeros);         // declared in analysis.h; defined in analysis.cpp

	NR=Realzeros.dimension();
	NC=Compzeros.dimension();
	Realzeros.sort();

	os+="\nTotal nunmber of roots: "+val::ToString(N);
	if (computerational) os+="\nNumber of rational roots found: " + val::ToString(Rationalzeros.length()) + " : \n";
	for (const auto &x : Rationalzeros ) os+=val::ToString(x) + "    ";
	os+="\nReal roots: "+val::ToString(NR) + " : \n";
	for (i=0;i<NR;i++) os+=val::ToString(Realzeros[i]) + "    ";

	os+="\nComplex roots: "+val::ToString(N-NR) + " , found: " + val::ToString(2*NC)+ " : \n";
	for (i=0;i<NC;i++) os+=val::ToString(Compzeros[i].real()) +"+-"+val::ToString(Compzeros[i].imaginary())+"i      ";
	
	if (f.degree() > 20 ) {
		std::cout << os << std::endl;
		return 0;
	}
	
	// Extrema:
	int NHP,NTP;
	val::vector<double> Maxima,Minima;

	N=val::extreme(f,Maxima,Minima);
	//std::cout<<os<<std::endl;
	if (N) {
		NHP=Maxima.dimension();
		NTP=Minima.dimension();
		os+="\n\nNumber of local extrema: " + val::ToString(N);
		if (NHP) {
			Maxima.sort();
			os+="\n\n Number of maxima: " + val::ToString(NHP) + "\n";
			for (i=0;i<NHP;i++) os+=" ( " + val::ToString(Maxima[i]) + " | " + val::ToString(f(Maxima[i])) + " )   ";
		}
		if (NTP) {
			Minima.sort();
			os+="\n\n Number of minima: " + val::ToString(NTP) + "\n";
			for (i=0;i<NTP;i++) os+=" ( " + val::ToString(Minima[i]) + " | " + val::ToString(f(Minima[i])) + " )   ";
		}
	}
	else {
		os+="\n\nNo local extrema.";
	}


	// Inflection Points:

	//val::pol<double> h(f);
	//h.derive();
	N=val::extreme(f.derive(),Maxima,Minima);

	if (N) {
		NHP=Maxima.dimension();
		NTP=Minima.dimension();
		os+="\n\nNumber of inflection points: " + val::ToString(N);
		if (NHP) {
			Maxima.sort();
			os+="\n\n Number of left-right: " + val::ToString(NHP) + "\n";
			for (i=0;i<NHP;i++) os+=" ( " + val::ToString(Maxima[i]) + " | " + val::ToString(f(Maxima[i])) + " )   ";
		}
		if (NTP) {
			Minima.sort();
			os+="\n\n Number of right-left: " + val::ToString(NTP) + "\n";
			for (i=0;i<NTP;i++) os+=" ( " + val::ToString(Minima[i]) + " | " + val::ToString(f(Minima[i])) + " )   ";
		}
	}
	else os+="\n\nNo inflection points.";

	std::cout<<os<<std::endl;

	/*
	if (argnr==3) {
		std::ofstream file(argv[2],std::ios::out | std::ios::trunc);
		if (!file) {
			std::cout<<"\nCannot write in "<<std::string(argv[2])<<"!!\n\n";
			return 0;
		}
		file<<os;
		std::cout<<"\nOutput written in "<<std::string(argv[2])<<"!!\n\n";
	}
	*/

	return 0;
}
