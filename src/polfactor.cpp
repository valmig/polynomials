#include <polfactor.h>
#include <Glist.h>
#include <fstream>
#include <val_utils.h>
#include <val_basics.h>
#include <valfunction.h>
#include <pol_arithmetic.h>
#include <MyTime.h>



// Actual algorithm is in <polfactor.h> and its source <polfactor.cpp> 

void program_info()
{
	std::cout<<"\nCommand: polfactor \"polynomial\"/Input-File [char = 0]";
	std::cout<<"\n=======";
	std::cout<<"\nComputes the irreducible factors of a polynomial over Q";
	std::cout<<"\nif char == 0, or over Fp, p = char, otherwise.\n";
}




int isinfix(const std::string &s)
{
    int n=s.size(),i;

    for (i=0;i<n;i++)
        if (s[i]=='x') return 1;
    return 0;

}


val::pol<val::modq> toModqPolynom(const val::pol<val::rational> & f)
{
	val::pol<val::modq> g;
	val::integer p(val::modq::q);
	val::polIterator<val::rational> monom;

	for (monom=f;monom;monom++) {
		g.insert(val::modq(val::nominator(monom()) % p),monom.actualdegree());
	}
	return g;
}



int main(int argnr, char* argv[])
{
	std::string firstarg;
	if (argnr >= 2) firstarg = std::string(argv[1]);
	
	if (argnr == 1 || (argnr == 2 && (firstarg == "-h" || firstarg == "--help") )) {
        program_info();
        return 0;
    }
	else if (argnr > 3) {
		std::cout << "\n Too many arguments. Quit program!!\n";
		return 0;
	}
    
	int cfield = 0;
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
	
	if (argnr == 3) cfield = val::FromString<int>(std::string(argv[2]));
	if (cfield < 0) cfield = 0;
	
	

	if (cfield) {  // polynomial in Fp[X]
		val::d_array<int> mult;
		val::pol<val::modq> f,h;
		//val::pol<val::rational> g;
		val::d_array<val::pol<val::modq> > faktor;
		int i,l;

		
		val::modq::q = cfield;
		f = ::toModqPolynom(f_r);
		
		//std::cout << "\n f_r = \n" << f_r;
		//std::cout << "\n f = \n" << f;

		if (f.iszero()) return 0;
		std::cout<<"\nComputing factors..."<<std::endl;
		l=polfactor(f,faktor,mult);
		std::cout<<"\nPolynomial factorized!"<<std::endl;
		std::cout<<"factots in F"<<val::modq::q<<" of f =\n"<<f<<std::endl;
		h=val::pol<val::modq>(val::modq(1),0);  // h =1
		for (i=0;i<l;i++) {
			std::cout<<std::endl<<faktor[i]<<"\n = "<<val::PolToString(faktor[i]);
			std::cout<<"\nMultiplicity: "<<mult[i]<<std::endl;
			h*=val::power(faktor[i],mult[i]);
		}
		if (f!=f.leader()*h) std::cout<<"Something went wrong!!!\n";
	}
	else {  // polynomial in Q[X]
		val::pol<val::rational> f, h(val::rational(1),0), fh;
		val::ChronoClass Chrono;

		f = std::move(f_r);
		

		if (f.iszero()) return 0;

		std::cout<<"\nCompute factors of f...:"<<std::flush;

		val::d_array<val::pol<val::rational>> factors = val::polfactor(f,1);
		std::cout<<"\n\nIrreducible factors computed!\nTime in s: "<<Chrono();

		// Compute Multiplicity of factors:
		val::d_array<int> mult(factors.length());
		int e, i = 0;

		fh = f;
		for (const auto &g : factors) {
			e = 0;
			while ((fh%g).iszero()) {
				++e;
				fh /= g;
			}
			mult[i] = e;
			++i;
			h *= val::power(g,e);
		}
		
		f = val::toRationalPolynom(val::primitivpart(f));
		
		if (f.LC().signum() != h.LC().signum()) h *= val::rational(-1);

		if (f!=h)  {
			std::cout<<"\nSomething went wrong!!!\n";
			std::cout<<"\n f = "<<val::PolToString(f);
			std::cout<<"\n h = "<<val::PolToString(h);
		}


		std::cout<<"\nFactors of f:\n";
		for (i = 0; i < factors.length(); ++i) {
			std::cout<<std::endl<<factors[i]<<"\n = "<<val::PolToString(factors[i])<<"\nMultiplicity: "<<mult[i]<<std::endl;
		}
	}

	std::cout<<std::endl;
	return 0;
}
