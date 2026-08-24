#include "global_declarations.h"

//Forward declaration of important functions used throught the program

template<typename T>
T **Allocate_2D(T **&m, int t1, int t2);

template<typename T>
T ***Allocate_3D(T ***&a, int IMAX, int JMAX, int KMAX);

template<class T>
void Delete_2D(T **&m, int d1, int d2);

template<class T>
void Delete_3D(T ***&m, int d1, int d2, int d3);

template<typename T>
T MinModLim(T U0, T U1, T U2, T dx);

template<typename T>
T ReconValue(T U0, T U1, T U2, T U3, T dx);

template<typename T>
T Sign(T a);

template<typename T>
T MinMod(T Ujpp, T Ujp, T Uj, T Ujm);

template<typename T>
T Max2(T a1, T a2);

template<typename T>
T Min2(T a1, T a2);

template<typename T>
T Max3(T a1, T a2, T a3);

template<typename T>
T Min3(T a1, T a2, T a3);

template<typename T>
T Max4(T a1, T a2, T a3, T a4);

template<typename T>
T Min4(T a1, T a2, T a3, T a4);

template<typename T>
T Log_Mean(T al, T ar);

template <typename T>
std::string ToStringWithPrecision(const T a_value, const int n = 16);

template <typename T>
std::string get_str_between_two_str(const T s, const T start_delim,
                                    const T stop_delim);

//Function definitions
template<typename T>
T **Allocate_2D(T **&m, int t1, int t2) {
    m = new T *[t1];
    for (int i = 0; i < t1; ++i) {
        m[i] = new T[t2];
        for (int j = 0; j < t2; ++j)
            m[i][j] = 0.0;
    }
    return m;
}

template<typename T>

T ***Allocate_3D(T ***&a, int IMAX, int JMAX, int KMAX) {
    a = new T **[IMAX + 1];

    for (int i = 0; i < IMAX + 1; i++) {

        a[i] = new T *[JMAX + 1];

        for (int j = 0; j < JMAX + 1; j++) {

            a[i][j] = new T[KMAX + 1];

            for (int k = 0; k < KMAX + 1; k++) {

                a[i][j][k] = 0.0;
            }
        }
    }
    return a;
}

template<class T>
void Delete_2D(T **&m, int d1, int d2) {

    for (int i = 0; i < d1; ++i) {
        delete[] m[i];
    }
    delete[] m;
}

template<class T>
void Delete_3D(T ***&m, int d1, int d2, int d3) {

    for (int i = 0; i < d1; i++) {
        for (int j = 0; j < d2; j++) {
            delete[] m[i][j];
        }

        delete[] m[i];
    }

    delete[] m;
}

template<typename T>
T MinModLim(T U0, T U1, T U2, T dx) {
    double Ur, Ul;

    Ur = (U2 - U1) / dx;
    Ul = (U1 - U0) / dx;
    if (Ur * Ul <= 0.0) return 0.0;
    if (fabs(Ur) < fabs(Ul)) return Ur;
    return Ul;
}

template<typename T>
T ReconValue(T U0, T U1, T U2, T U3, T dx) {

    return (U2 - 0.5 * MinModLim(U1, U2, U3, dx) * dx) -
           (U1 + 0.5 * MinModLim(U0, U1, U2, dx) * dx); //For second order accuracy.
    //return (U2-U1); //For first order accuracy
}

template<typename T>
T Sign(T a) {

    if (a > 0.0) return 1.0;
    else if (a < 0.0) return -1.0;
    else return 0.0;
}

/*
template <typename T>
T MinMod(T Ujpp, T Ujp, T Uj, T Ujm){
	double r1,r2,del=0.000001;

	if (fabs(Ujp-Uj)<del) r1=Sign(Ujp-Uj)*(Uj-Ujm)/del;
	else r1=(Uj-Ujm)/(Ujp-Uj);

	if (fabs(Ujp-Uj)<del) r2=Sign(Ujp-Uj)*(Ujpp-Ujp)/del;
	else r2=(Ujpp-Ujp)/(Ujp-Uj);
	
	if (r1*r2<=0) return 0.0;
	else if ((fabs(r1)<fabs(r2)) and r1*r2<1.0) return r1;
	else if ((fabs(r2)<fabs(r1)) and r1*r2<1.0) return r2;
	else return 1.0;
}
*/


template<typename T>
T MinMod(T del_plus, T delu, T del_minus) {
    double r1, r2, del = 0.000001;

    if (fabs(delu) < del) r1 = Sign(delu) * (del_minus) / del;
    else r1 = (del_minus) / (delu);

    if (fabs(delu) < del) r2 = Sign(delu) * (del_plus) / del;
    else r2 = (del_plus) / (delu);

    if (r1 * r2 <= 0) return 0.0;
    else if ((fabs(r1) < fabs(r2)) and (fabs(r1) < 1.0) and r1 * r2 > 0.0) return r1;
    else if ((fabs(r2) < fabs(r1)) and (fabs(r2) < 1.0) and r1 * r2 > 0.0) return r2;
    else return 1.0;
}

template<typename T>
T Max3(T a1, T a2, T a3) {

    if (a1 >= a2 && a1 >= a3) return a1;
    if (a2 >= a3 && a2 >= a1) return a2;
    else return a3;

}

template<typename T>
T Max2(T a1, T a2) {

    if (a1 >= a2) return a1;
    else return a2;

}

template<typename T>
T Min2(T a1, T a2) {

    if (a1 >= a2) return a2;
    else return a1;

}


template<typename T>
T Min3(T a1, T a2, T a3) {

    if (a1 <= a2 && a1 <= a3) return a1;
    if (a2 <= a3 && a2 <= a1) return a2;
    else return a3;

}


template<typename T>
T Max4(T a1, T a2, T a3, T a4) {

    if (a1 >= Max3(a2, a3, a4)) return a1;
    else return Max3(a2, a3, a4);

}

template<typename T>
T Min4(T a1, T a2, T a3, T a4) {

    if (a1 <= Min3(a2, a3, a4)) return a1;
    else return Min3(a2, a3, a4);

}


template<typename T>
T Log_Mean(T al, T ar) {
    double xi, ln_xi, f, u, eps = 0.01, F;

    xi = al / ar;
    f = (xi - 1.0) / (xi + 1.0);
    u = f * f;

    if (u < eps) {
        /*series approximation of 2*f/ln(xi), accurate for near-equal states*/
        F = 1.0 + u / 3.0 + u * u / 5.0 + u * u * u / 7.0;
        return (al + ar) / 2.0 * F;
    } else {
        /*exact logarithmic mean; the previous expression inverted f/ln(xi)
          and lost the sign, returning negative values whenever al/ar > ~1.1
          or < ~0.9 outside the series range*/
        ln_xi = std::log(xi);
        return (al - ar) / ln_xi;
    }
}

template<typename T>
T Minmod_Flux(T del_p, T del_m) {

    double r, eps = 1e-6;

    r = 0.0;

    r = del_p / (del_m + 1e-16);

    if (fabs(del_p)<eps) {

        r = Sign(del_p)*del_m/eps;

    }

    return std::max(0.0, std::min(1.0, r));
    //return r;
}

template<typename T>
T VanAlbada_Limiter(T del_plus, T del_minus) {

    const double omega = 1e-12;
    double r;

    //del_minus = Sign(del_minus)*(fabs(del_minus)+omega);

    r = del_plus / (del_minus + 1e-15);

    return ((r * r + r) / (r * r + 1.0));
}

template<typename T>
T Venki_Limiter(T del_plus, T del_minus, T dx) {

    double const k = 0.1;
    double eps;

    eps = pow(k * dx, 3);

    return ((del_plus * del_plus + eps) * del_minus + (del_minus * del_minus + eps) * del_plus) /
           (del_plus * del_plus + del_minus * del_minus + 2.0 * eps + 1e-10);
}

template <typename T>
std::string ToStringWithPrecision(const T a_value, const int n){
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

template <typename T>
std::string get_str_between_two_str(const T s, const T start_delim,
                                    const T stop_delim){
    unsigned first_delim_pos = s.find(start_delim);
    unsigned end_pos_of_first_delim = first_delim_pos + start_delim.length();
    /*search the closing delimiter AFTER the opening one; searching from the
      start returned the first occurrence again and yielded everything past
      the opening delimiter instead of the text between them*/
    unsigned last_delim_pos = s.find(stop_delim, end_pos_of_first_delim);

    if (last_delim_pos == std::string::npos)
        return s.substr(end_pos_of_first_delim);

    return s.substr(end_pos_of_first_delim,
                    last_delim_pos - end_pos_of_first_delim);
}