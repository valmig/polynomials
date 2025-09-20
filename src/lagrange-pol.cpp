#include <iostream>
#include <val_utils.h>
#include <d_array.h>
#include <rational.h>
#include <pol.h>
#include <pol_arithmetic.h>


void program_info()
{
    std::cout<<"\nCommand: lagrange-pol \"nodes\"";
    std::cout<<"\n=======";
    std::cout<<"\nComputes the lagrange basis polynomials wrt the nodes.\n";
}


val::d_array<val::pol<val::rational>> lagrangepols(const val::d_array<val::rational> &values)
{
    val::d_array<val::pol<val::rational>> P;

    if (values.isempty()) return  P;


    val::d_array<val::rational> rvalues = values;
    int i, j, m = rvalues.length();
    val::pol<val::rational> X(1,1), one(1,0), l, p;

    rvalues.sort();

    for (i = 0; i < m; ++i) {
        if (i != 0 && rvalues[i-1] == rvalues[i]) continue;
        p.del();
        p = one;
        for (j = 0; j < m ; ++j) {
            if (j != 0 && rvalues[j-1] == rvalues[j]) continue;
            if (i == j) continue;
            l = X; l.insert(-rvalues[j],0); l /= rvalues[i]-rvalues[j];
            p *= l;
        }
        P.push_back(std::move(p));
    }

    return P;
}



int main(int argnr, char *argv[])
{
    std::string firstarg;
    if (argnr == 2) firstarg = std::string(argv[1]);

    if (argnr == 1 || (argnr == 2 && (firstarg == "-h" || firstarg == "--help") )) {
        program_info();
        return 0;
    }

    val::d_array<std::string> svalues = val::getwordsfromstring<char, val::d_array, val::d_array>(firstarg,val::d_array<char>{' '});
    val::d_array<val::rational> rvalues;

    rvalues.reserve(svalues.length());

    for (const auto &s : svalues) {
        rvalues.push_back(val::FromString<val::rational>(s));
    }
    val::d_array<val::pol<val::rational>> P = lagrangepols(rvalues);

    std::cout << "\n Lagrange Polynomials:\n";
    for (const auto &p : P) {
        std::cout << val::PolToString(p) << std::endl;
    }

    std::cout << std::endl;
    return 0;
}
