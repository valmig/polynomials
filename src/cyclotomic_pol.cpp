#include <rational.h>
#include <numbers.h>
#include <pol.h>
#include <pol_arithmetic.h>



void program_info()
{
	std::cout<<"\nCommand: cyclotomic_pol n";
	std::cout<<"\n=======";
	std::cout<<"\nComputes the cyclotomic polynomial of degree n in Z[x]\n";
}

// p must be a prime number
val::pol<val::integer> cyclotomic_prime(int p)
{
	val::pol<val::integer> f;
	val::integer one(1);
	
	for (int i = 0; i < p; ++i) {
		f.insert(one,i);
	}
	
	return f;
}

val::pol<val::integer> cyclotomic_pol(int n)
{
	if ( n<= 0) return val::pol<val::integer>();
	
	val::integer one(1);
	val::pol<val::integer> f({ {one,1},{-one,0} }), X;
	int m = 1, n1 = n, p;
	
	
	if (n == 1) return f;
	
	for (p = 2; p <= n1; p = val::nextprime(p)) {
		if (n1 % p) continue;
		m *= p;
		while ( !(n1 % p) ) {
			n1/=p;
		}
		X = val::pol<val::integer>(one,p);
		f = f(X)/f;
	}
	n /= m;
	if (n!=1) {
		X = val::pol<val::integer>(one,n);
		f = f(X);
	}
	
	return f;	
}

void test()
{
	int n;
	val::pol<val::integer> f;
	
	do {
		std::cout << "\ninput n : ";
		std::cin >> n;
		f = cyclotomic_pol(n);
		std::cout << " f = " << val::PolToString(f) << std::endl;
	}
	while (n);
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
		return 1;
	}
	
	val::pol<val::integer> f = cyclotomic_pol(val::FromString<int>(firstarg));
	
	std::cout << "\n" << f << "\n = " << val::PolToString(f);
	
	std::cout << std::endl;
	return 0;
}
