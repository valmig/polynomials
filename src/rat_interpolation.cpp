#include <iostream>
#include <rational.h>
#include <pol.h>
#include <string>
#include <val_utils.h>
#include <d_array.h>
#include <LA.h>
#include <pol_arithmetic.h>


void program_info()
{
 std::cout<<"\nCommand: rat_interpolation \"pairs of numbers\"";
 std::cout<<"\n=======";
 std::cout<<"\nComputes the interpolation-polynomial over Q wrt the given pairs!\n";
}

int isinarray(const val::rational& value,const val::d_array<val::rational> & x)
{
    for (const auto& v : x)
        if (v==value) return 1;
    return 0;
}


val::matrix<val::rational> set_les(const std::string &s)
{
    int n=s.length(),i,j;
    val::matrix<val::rational> A;
    if (n==0) return A;
    val::d_array<val::rational> x,y;
    val::d_array<val::rational>::set_plus_cap(10);
    std::string svalue="";
    val::rational value;

    j=0;
    for (i=0;i<n;++i) {
        if (s[i]==' ' || s[i]== '\n') {
            j%=2;
            value=val::FromString<val::rational>(svalue);
            svalue="";
            if (!j) {
                if (!isinarray(value,x)) {
                    x.push_back(value);
                    ++j;
                }
            }
            else {
                y.push_back(value);
                ++j;
            }
        }
        else svalue+=s[i];
    }

    value=val::FromString<val::rational>(svalue);
    if (!j) {
        if (!isinarray(value,x)) {
            x.push_back(value);
            ++j;
        }
    }
    else {
        y.push_back(value);
        ++j;
    }

    n= val::Min(x.length(),y.length());
    if (n==0) return A;
    A = val::matrix<val::rational>(n,n+1);
    for (i=0;i<n;++i)
        for (j=0;j<n;++j) A(i,j) = val::power(x[i],j);
    for (i=0;i<n;++i) A(i,n) = y[i];
    return A;
}

val::matrix<val::rational> set_les2(const std::string &s)
{
    int n=s.length(),i=0,j,k=0,n1=0,n2=0,l;
    val::matrix<val::rational> A;
    if (n==0) return A;
    val::d_array<val::d_array<val::rational>> x(3),y(3);
    val::d_array<val::rational>::set_plus_cap(10);
    val::d_array<int> dim(0,3);
    std::string svalue="";
    val::rational value;

    while (s[i]==' ' || s[i]=='\n' || s[i] == ';') ++i;
    j=0;
    for (;i<n;++i) {
        if (s[i]==' ' || s[i]== '\n' || s[i]==';') {
            j%=2;
            value=val::FromString<val::rational>(svalue);
            svalue="";
            if (!j) {
                //if (!isinarray(value,x[k])) {
                if (!val::isinContainer(value,x[k])) {
                    x[k].push_back(value);
                    ++j;
                }
            }
            else {
                y[k].push_back(value);
                ++j;
            }
            if (s[i]==';' && k<2) k++;
        }
        else svalue+=s[i];
    }

    value=val::FromString<val::rational>(svalue);
    if (!j) {
        if (!isinarray(value,x[k])) {
            x[k].push_back(value);
            ++j;
        }
    }
    else {
        y[k].push_back(value);
        ++j;
    }

    dim[1]=n1= val::Min(x[1].length(),y[1].length());
    dim[2]=n2= val::Min(x[2].length(),y[2].length());
    dim[0]=n= val::Min(x[0].length(),y[0].length());
    n+=n1+n2;
    if (n2) n = val::Max(n,n2+2);
    if (n1) n = val::Max(n,n1+1);
    if (n==0) return A;
    /*
    for (const auto & v : x) {
        for (const auto &wert : v) std::cout<<wert<<"  ";
        std::cout<<std::endl;
    }
    */
    //std::cout<<n<<"  "<<dim[0]<<"  "<<dim[1]<<"  "<<dim[2]<<std::endl;
    A = val::matrix<val::rational>(n,n+1);
    for (k=0,l=0;k<3;++k) {
        for (i=0;i<dim[k];++i,++l) {
            for (j=k;j<n;++j) {
                A(l,j) = val::power(x[k][i],j-k);
                if (k>=1) A(l,j) *= val::rational(j);
                if (k==2) A(l,j) *= val::rational(j-1);
            }
            A(l,n) = y[k][i];
        }
    }
    return A;
}


val::pol<val::rational> interpolation(const std::string &s)
{
    val::matrix<val::rational> X,A=set_les2(s);
    val::pol<val::rational> f;
    val::rational det;
    int dim;

    if (A.numberofcolumns()<2) return f;
    //std::cout<<"\nMatrix: \n"<<A<<std::endl;
    dim=val::les(A,X,det);
    if (dim==0) return f;
    for (int i=0;i<X.numberofcolumns();++i) f.insert(X(0,i),i);
    return f;
}


void test_interpolation()
{
    std::string s = "-2 0 2 0;-2 -8 2 -8;-2 0";
    std::cout<<"\n A = \n"<<set_les2(s);
    std::cout<<"\n f = \n"<<interpolation(s);
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
    
    val::pol<val::rational> f = interpolation(firstarg);

    std::cout << "\nInterpolation-polynomial:\n" << f << "\n = " << val::PolToString(f) << std::endl;
    return 0;
}
