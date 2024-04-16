#include <polfactor.h>
#include <valfunction.h>
#include <pol_arithmetic.h>
#include <LA.h>
#include <fraction.h>
#include <error.h>



void program_info()
{
 std::cout<<"\nCommand: partial-fraction \"rational function\" ";
 std::cout<<"\n=======";
 std::cout<<"\nTries to compute the partial fraction decomposition of the given rational function.\n";
}




int sortfactors(val::d_array<val::pol<val::rational>> &factors, int &r, int &s)
{
    using namespace val;
    r = s = 0;
    for (const auto & f : factors) {
        if (f.degree() == 1) ++r;
        else if (f.degree() == 2) ++s;
        else return 0;
    }
    if (s == 0 || r == 0) return 1;
    int n = factors.length(),j;
    d_array<pol<rational>> hfactors(n);

    j=0;
    for (auto &f : factors) {
        if (f.degree() == 1) {
            hfactors[j] = std::move(f);
            ++j;
        }
    }
    j=0;
    for (auto &f : factors) {
        if (f.degree() == 2) {
            hfactors[j+r] = std::move(f);
            ++j;
        }
    }
    factors = std::move(hfactors);

    return 1;
}


int partialfraction(const val::valfunction& f, val::rational &cont, val::pol<val::rational> &fp,val::d_array<val::pol<val::rational>> &numpol,
                    val::d_array<val::pol<val::rational>> &denompol, val::d_array<int> &denumexpo, int comment = 0)
{
    using namespace val;
    numpol.del(); denompol.del(); denumexpo.del(); fp.del();

    if (!f.isrationalfunction()) return 0;
    rationalfunction fR = f.getrationalfunction();
    pol<rational> p = fR.nominator(), q = fR.denominator(), h, hdiv;

    if (deg(p) >= deg(q)) {
        pol<rational> rpol;
        divrem(p,q,fp,rpol);
        if (p != q*fp + rpol) {
            Error::error("\n Error by computation of division with remainder!");
        }
        p = std::move(rpol);
    }

    cont = content(p)/content(q);
    p = toRationalPolynom(primitivpart(p));
    q = toRationalPolynom(primitivpart(q));

    d_array<pol<rational>> qfactors = polfactor(q);
    int i, j, k, n = qfactors.length(), r, s, m = 0, nqvar = 0, kb;
    d_array<int> qexponents(0,n);
    d_array<int> islinear(0,n);

    // set linear factors first, r is number of linear factors, s of quadratic factors.
    if (!sortfactors(qfactors,r,s)) return 0;

    h = q;
    for (i = 0; i < n; ++i) {
        while ((h%qfactors[i]).iszero()) {
           ++qexponents[i];
           h /= qfactors[i];
        }
        if (i<r) m+=qexponents[i];
        else {
            m+=2*qexponents[i];
            nqvar += qexponents[i];
        }
    }
    if (comment) std::cout<<"\n After factorization: m = "<<m<<", h = "<<PolToString(h)<< ", cont = "<<ToString(cont);
    /*
    std::cout<<std::endl;
    for (const auto& v : qfactors) std::cout<<PolToString(v)<<std::endl;
    */
    d_array<d_array<pol<rational>>> linq(r), quadq(s);
    d_array<d_array<int>> a_var(r), alpha_var(s), beta_var(s);
    for (i = 0; i < r; ++i) {
        linq[i] = d_array<pol<rational>>(qexponents[i]);
        a_var[i] = d_array<int>(qexponents[i]);
    }
    for (i = 0; i < s; ++i) {
        quadq[i] = d_array<pol<rational>>(qexponents[i+r]);
        alpha_var[i] = d_array<int>(qexponents[i+r]);
        beta_var[i] = d_array<int>(qexponents[i+r]);
    }
    k = 0;

    for (i = 0; i < r; ++i) {
        hdiv = pol<rational>(rational(1),0);
        for (j = 0; j < qexponents[i]; ++j) {
            hdiv *= qfactors[i];
            linq[i][j] = q/hdiv;
            a_var[i][j] = k;
            ++k;
        }
    }
    for (i = 0; i < s; ++i) {
        hdiv = pol<rational>(rational(1),0);
        for (j = 0; j < qexponents[i+r]; ++j) {
            hdiv *= qfactors[i+r];
            quadq[i][j] = q/hdiv;
            alpha_var[i][j] = k;
            beta_var[i][j] = k + nqvar;
            ++k;
        }
    }
    //if (k+ nqvar != m) std::cout<<"\n k+nqvar != m! k + nqvar = "<<k+nqvar;

    // Set LES:
    matrix<rational> A(m,m+1);

    for (int d = 0; d < m; ++d) {
        A(d,m) = p[d];  // last column:
        // a_var
        for (i = 0; i < r; ++i) {
            for (j = 0; j < qexponents[i]; ++j) {
                k = a_var[i][j];
                A(d,k) = linq[i][j][d];
            }
        }
        // alpha_var, beta_var
        for (i = 0; i < s; ++i) {
            for (j = 0; j < qexponents[i+r]; ++j) {
                k = alpha_var[i][j];
                if (d > 0) A(d,k) = quadq[i][j][d-1];
                k = beta_var[i][j];
                A(d,k) = quadq[i][j][d];
            }
        }
    }

    if (comment) {
        for (i = 0; i < m; ++i) {
            std::cout<<std::endl;
            for (j = 0; j <= m; ++j) std::cout<<A(i,j)<<"  ";
        }
    }
    matrix<rational> X;
    rational det;
    les(A,X,det);
    if (X.numberofrows() != 1) return 0;

    numpol.resize(m - nqvar); denompol.resize(m - nqvar); denumexpo.resize(m -nqvar);


    for (i = 0; i < r; ++i) {
        for (j = 0; j < qexponents[i]; ++j) {
            k = a_var[i][j];
            numpol[k] = pol<rational>(X(0,k),0);
            denompol[k] = qfactors[i];
            denumexpo[k] = j+1;
        }
    }
    for (i = 0; i < s; ++i) {
        for (j = 0; j < qexponents[i+r]; ++j) {
            k = alpha_var[i][j];
            kb = beta_var[i][j];
            numpol[k] = pol<rational>({GPair<rational,int>(X(0,k),1),GPair<rational,int>(X(0,kb),0) });
            denompol[k] = qfactors[i+r];
            denumexpo[k] =  j+1;
        }
    }

    //delete 0-entries:
    n = 0 , m = numpol.length();
    for (i = 0; i < m; ++i) if (!numpol[i].iszero()) ++n;
    d_array<pol<rational>> hnumpol(n), hdenompol(n);
    d_array<int> hdenumexpo(n);

    j = 0;
    for (i = 0; i < m; ++i) {
        if (!numpol[i].iszero()) {
            hnumpol[j] = std::move(numpol[i]);
            hdenompol[j] = std::move(denompol[i]);
            hdenumexpo[j] = denumexpo[i];
            ++j;
        }
    }

    numpol = std::move(hnumpol);
    denompol = std::move(hdenompol);
    denumexpo = std::move(hdenumexpo);

    return 1;
}



val::d_array<val::valfunction> partialfraction(const val::valfunction& f)
{
    using namespace val;
    d_array<valfunction> F;

    if (!f.isrationalfunction()) return F;
    rationalfunction fR = f.getrationalfunction();
    pol<rational> p = fR.nominator(), q = fR.denominator(), fp, h, hdiv;

    if (deg(p) >= deg(q)) {
        pol<rational> rpol;
        divrem(p,q,fp,rpol);
        if (p != q*fp + rpol) {
            std::cout<<"\n Error by computation of division with remainder";
            exit(-1);
        }
        p = std::move(rpol);
        //std::cout<<"\n After div with rem: f(x) = " + PolToString(fp) + " + (" + PolToString(p) + ")/(" + PolToString(q) + ")";
    }

    rational cont = content(p)/content(q);
    p = toRationalPolynom(primitivpart(p));
    q = toRationalPolynom(primitivpart(q));

    d_array<pol<rational>> qfactors = polfactor(q);
    int i, j, k, n = qfactors.length(), r, s, m = 0, nqvar = 0, kb;
    d_array<int> qexponents(0,n);
    d_array<int> islinear(0,n);

    // set linear factors first, r is number of linear factors, s of quadratic factors.
    if (!sortfactors(qfactors,r,s)) return F;
    std::cout<<"\n r = "<<r<<" , s = "<<s;

    h = q;
    for (i = 0; i < n; ++i) {
        while ((h%qfactors[i]).iszero()) {
           ++qexponents[i];
           h /= qfactors[i];
        }
        if (i<r) m+=qexponents[i];
        else {
            m+=2*qexponents[i];
            nqvar += qexponents[i];
        }
    }
    std::cout<<"\n After factorization: m = "<<m<<", h = "<<PolToString(h);
    /*
    std::cout<<std::endl;
    for (const auto& v : qfactors) std::cout<<PolToString(v)<<std::endl;
    */
    d_array<d_array<pol<rational>>> linq(r), quadq(s);
    d_array<d_array<int>> a_var(r), alpha_var(s), beta_var(s);
    for (i = 0; i < r; ++i) {
        linq[i] = d_array<pol<rational>>(qexponents[i]);
        a_var[i] = d_array<int>(qexponents[i]);
    }
    for (i = 0; i < s; ++i) {
        quadq[i] = d_array<pol<rational>>(qexponents[i+r]);
        alpha_var[i] = d_array<int>(qexponents[i+r]);
        beta_var[i] = d_array<int>(qexponents[i+r]);
    }
    k = 0;

    for (i = 0; i < r; ++i) {
        //std::cout<<std::endl;
        hdiv = pol<rational>(rational(1),0);
        for (j = 0; j < qexponents[i]; ++j) {
            hdiv *= qfactors[i];
            linq[i][j] = q/hdiv;
            a_var[i][j] = k;
            ++k;
            //std::cout<<a_var[i][j]<<"  ";
        }
    }
    for (i = 0; i < s; ++i) {
        //std::cout<<std::endl;
        hdiv = pol<rational>(rational(1),0);
        for (j = 0; j < qexponents[i+r]; ++j) {
            hdiv *= qfactors[i+r];
            quadq[i][j] = q/hdiv;
            alpha_var[i][j] = k;
            beta_var[i][j] = k + nqvar;
            //std::cout<<alpha_var[i][j]<<" , "<<beta_var[i][j]<<" , ";
            ++k;
        }
    }
    if (k+ nqvar != m) std::cout<<"\n k+nqvar != m! k + nqvar = "<<k+nqvar;
    // Set a_var, alpha_va, beta_var;


    // Set LES:

    matrix<rational> A(m,m+1);

    for (int d = 0; d < m; ++d) {
        A(d,m) = p[d];  // last column:
        // a_var
        for (i = 0; i < r; ++i) {
            for (j = 0; j < qexponents[i]; ++j) {
                k = a_var[i][j];
                A(d,k) = linq[i][j][d];
            }
        }
        // alpha_var, beta_var
        for (i = 0; i < s; ++i) {
            for (j = 0; j < qexponents[i+r]; ++j) {
                k = alpha_var[i][j];
                if (d > 0) A(d,k) = quadq[i][j][d-1];
                k = beta_var[i][j];
                A(d,k) = quadq[i][j][d];
            }
        }
    }

    for (i = 0; i < m; ++i) {
        std::cout<<std::endl;
        for (j = 0; j <= m; ++j) std::cout<<A(i,j)<<"  ";
    }
    matrix<rational> X;
    rational det;
    les(A,X,det);
    if (X.numberofrows() != 1) return F;
    //std::cout<<std::endl;
    //for (i = 0; i < m; ++i) std::cout<<X(0,i)<<"  ";
    if (fp.iszero()) F = d_array<valfunction>(m - nqvar);
    else {
        F = d_array<valfunction>(m - nqvar +1);
        F[m - nqvar] = valfunction(PolToString(fp));
    }
    d_array<int> denomexp(m-nqvar);

    std::string scont = ToString(cont),sf;

    for (i = 0; i < r; ++i) {
        for (j = 0; j < qexponents[i]; ++j) {
            k = a_var[i][j];
            denomexp[k] = j + 1;
            sf = "(" + scont + ")*" + ToString(X(0,k)) + "/ (" + PolToString(qfactors[i]) + ")^" + ToString(denomexp[k]); //ToString(qexponents[i]);
            //std::cout<<std::endl<<sf;
            //F[k] = valfunction("(" + scont + ")*" + ToString(X(0,k)) + "/ (" + PolToString(qfactors[i]) + ")^" + ToString(qexponents[i]));
            F[k] = valfunction(sf);
        }
    }
    //std::cout<<"\n OK!"<<std::endl;
    for (i = 0; i < s; ++i) {
        for (j = 0; j < qexponents[i+r]; ++j) {
            k = alpha_var[i][j];
            kb = beta_var[i][j];
            denomexp[k] = j + 1;
            sf = "(" + scont + ")*(" + ToString(X(0,k)) + "*x + " + ToString(X(0,kb)) + ")/ (" + PolToString(qfactors[i+r]) + ")^" + ToString(denomexp[k]); //ToString(qexponents[i+r]);
            //std::cout<<std::endl<<sf;
            F[k] = valfunction(sf);
        }
    }

    return F;
}



void test_partialfraction()
{
    //val::valfunction f("(x -7)/((x+1)(x-2)^2(x^2+1)^2)");
    //val::valfunction f("(1-x)/((x+1)(x^2 +1))");
    //val::valfunction f("(5x^2 - x +8)/((x+2)(x^2 -2x+2))");
    //val::valfunction f("(4x-1)/(x^2 + x -2)");
    //val::valfunction f("(2x-3)/((x-1)^2(x+1))");
    val::valfunction f("(1-x)/(x^2-1)");
    auto F = partialfraction(f);

    std::cout<<"\n Partial fractions: \n";
    for (const auto &v : F) {
        std::cout<<v.getinfixnotation()<<" , ";
    }
}


void outputresult(const std::string &sf)
{
    using namespace val;
    valfunction f(sf);
    integer i_one(1);
    rational cont, one(1), lc;
    pol<rational> fp;
    d_array<int> denomexp;
    d_array<pol<rational>> numpol, denompol;
    if (!partialfraction(f,cont,fp,numpol,denompol,denomexp,1)) {
        std::cout<<"\nPartial fraction decomposition failed!";
        return;
    }
    int n = numpol.length(), contpar = 0, i;
    std::string os;
    if (!fp.iszero()) {
        os += PolToString(fp);
        if (n > 0) os += " + ";
    }
    if (!n) return;
    if (cont == -one) {
        os += "-";
        if (n>1) {
            os += "(";
			contpar = 1;
        }
    }
    else if (cont != one) {
        os += ToString(cont) + " * ";
        if (n>1) {
            os +="(";
            contpar = 1;
        }
    }

    for (i = 0; i < n; ++i) {
        if (numpol[i].degree() == 0) {
            lc = numpol[i].LC();
            if (lc.signum() < 0 || lc.denominator() != i_one) os += "(" + ToString(lc) + ")";
            else os += ToString(lc);
        }
        else if (numpol[i].length() > 1) {
            os += "(" + PolToString(numpol[i]) + ")";
        }
        else {
            if (numpol[i][0] != one) os += PolToString(numpol[i]);
            else os += "1";
        }
        os+="/";
        if (denompol[i].length() > 1) os += "(" + PolToString(denompol[i]) + ")";
        else os += PolToString(denompol[i]);
        if (denomexp[i] > 1) os += "^" + ToString(denomexp[i]);
        if (i != n-1) os += " + ";
    }
    if (contpar) os += ")";
    std::cout<<"\nPartial fraction decomposition: \n"<<os;
}



int main(int argnr,char* argv[])
{
    //test_partialfraction();

    if (argnr == 1 || (argnr > 1 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help"))) {
        program_info();
        return 0;
    }
    outputresult(std::string(argv[1]));

    std::cout<<std::endl;
    return 0;
}
