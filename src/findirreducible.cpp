#include <rand.h>
#include <modq.h>
#include <integer.h>
#include <pol.h>
#include <string>
#include <val_utils.h>
#include <pol_arithmetic.h>




void program_info()
{
 std::cout<<"\nCommand: findirreducible p n";
 std::cout<<"\n=======";
 std::cout<<"\nCreates a random irreducible polynomial of degree n in Fp\n";
}


// Creates random polynomial of degree n
val::pol<val::modq> randpoly(int n)
{

 int i;
 val::pol<val::modq> h;

 for (i=0;i<n;++i) h.insert(val::modq(val::random(val::modq::q)),i);
 h.insert(val::modq(1),n);           // Leitcoeffz = 1
 return h;
}


// Computes (a^n mod basis)  non-recursive.
template <class T>
val::pol<T> powermod(val::pol<T> a,val::integer n,const val::pol<T>& basis)
{
 val::pol<T> x(T(1),0);
 val::integer zwei(2);

 while (n!=0) {
      if (!n.iseven()) {
         x*=a;
         if (x.degree()>=basis.degree()) x%=basis;
       }
       a*=a;
       if (a.degree()>=basis.degree()) a%=basis;
       n/=zwei;
 }
 return x;
}



val::pol<val::modq> findirreducible(int n)
{
 int i,m,fertig;
 val::pol<val::modq> f,g,x(val::modq(1),1); // x= x
 val::integer q(val::modq::q);

 if (n==1) return randpoly(1);
 m=n/2;

 fertig=0;
 while (!fertig) {
       //randomize();
       f=randpoly(n);
       if (f[0]==val::modq(0)) f.insert(val::modq(1),0);
       g=powermod(x,q,f);
       for (i=1;i<=m;i++)
	   if (val::gcd(g-x,f)!=1) {
	      fertig=0;
	      break;
	   }
	   else {
	      fertig=1;
	      g=powermod(g,q,f);
	   }
 }
 return f;
}


int main(int argnr,char* argv[])
{
	std::string firstarg;
	if (argnr >= 2) firstarg = std::string(argv[1]);
	
	if (argnr == 1 || (argnr == 2 && (firstarg == "-h" || firstarg == "--help") )) {
        program_info();
        return 0;
    }
    if (argnr != 3) {
		std::cout << "\nProgram requires 2 arguments!" << std::endl;
		return 1;
	}
    
    val::initialize_random();
    val::modq::q = val::FromString<int>(firstarg);
    int n = val::FromString<int>(std::string(argv[2]));
    if (n<=0) return 0;
    
    val::pol<val::modq> f = findirreducible(n);

    std::cout<<"\n Irreducible polynomial of degree "<<n<<" in F"<<val::modq::q<<":\n" << f;
    std::cout << "\n = " << val::PolToString(f) << std::endl;


    return 0;
}
