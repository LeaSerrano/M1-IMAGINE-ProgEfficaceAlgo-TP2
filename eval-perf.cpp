#include "eval-perf.h"
#include <iostream>
#include <x86intrin.h>
#include <chrono>
#include <vector>
#include <math.h>
#include <cmath>

#define C1 0.2f
#define C2 0.3f

void Ex4_test() {
    long int n = 0;
    for (long int i = 0; i < 100000; i++) {
        n *= i;
    }
}

std::vector<int> Ex5_sommePrefixe(std::vector<int> tab, std::vector<int> tabFinal) {

    for (int i = 0; i < tab.size(); i++) {
        if (i == 0) {
            tabFinal[i] = tab[i];
        }
        else {
            tabFinal[i] = tab[i] + tab[i-1];
        }
    }

    //Essai d'optimisation non fonctionnelle
   /*for (int i = 0; i < milieuTab; i++) {
        if (i > 0) { 
            tabP1[i] = tabP1[i-1] + tab[i];
        }
        else {
            tabP1[i] = tab[i];
        }

        //std::cout << tabP1[i] << " " << i << std::endl;
    }

    for (int j = milieuTab; j < tab.size(); j++) {
        if (j > milieuTab) {
            tabP2[j] = tabP2[j-1] + tab[j];
        }
        else {
            tabP2[j] = tab[j];
        }

        //std::cout << tabP2[j] << " " << j << std::endl;
    }

    for (int k = 0; k < tab.size(); k++) {
        if (k < milieuTab) {
            tabFinal[k] = tabP1[k];
        }
        else {
            tabFinal[k] = tabFinal[milieuTab-1] + tabP2[k];
        }

        std::cout << tabFinal[k] << " " << k << std::endl;
    }*/

    return tabFinal;
}

float Ex6_puissances_alpha (float resultat, std::vector<float> P, int x) {
    resultat = 0;
    for (int i = 0; i < P.size(); i++) {
        resultat += P[i]*(pow(x, i));
    }
    return resultat;
}

float Ex6_Horner (float resultat, std::vector<float> P, int x) {
    for (int i = P.size()-1; i > 0 ; i--) {
        if (i == (P.size()-1)) {
            resultat = P[i];
        }
        else {
            resultat = resultat*x + P[i];
        }
    }
    return resultat;
}

std::vector<float> Ex7_reduce_mult(std::vector<std::vector<float>> V , std::vector<float> res) {
    for (int j = 0; j < V[0].size(); j++) {
        for (int i = 0; i < V.size(); i++) {
            res[j] = res[j] * V[i][j];
        }
    }
    return res;
}//CPI = 12045

std::vector<float> Ex7_reduce_plus(std::vector<std::vector<float>> V , std::vector<float> res) {
    for (int j = 0; j < V[0].size(); j++) {
        for (int i = 0; i < V.size(); i++) {
            res[j] = res[j] + V[i][j];
        }
    }
    return res;
}//CPI = 11800

//Ex8
void slowperformance1(std::vector<float> x, std::vector<float> y, std::vector<float> z, int n) {
    for(int i = 0 ; i < (n-2) ; i ++) {
        x[i] = x[i]/M_SQRT2 + y[i] * C1;
        x[i+1] += z[(i%4) * 10] * C2;
        x[i+2] += sin((2 * M_PI * i)/3 ) * y[i+2]; 
        
    }
}//nbCycles ~= 2060289

void slowperformance2(std::vector<float> x, std::vector<float> y, std::vector<float> z, int n) {
    x[0] = x[0]/M_SQRT2 + y[0] * C1;
    x[1] = x[1] + (x[1]/M_SQRT2 + y[1] * C1);

    for(int i = 2 ; i < (n-2) ; i ++) {
        x[i] = x[i]/M_SQRT2 + y[i] * C1 + z[((i-1)%4) * 10] * C2 + sin((2 * M_PI * (i-2))/3 ) * y[i];
    }

}//nbCycles ~= 1913491


float Ex9_puissances_alpha (float resultat, std::vector<float> P, int x, int i, int sizeP) {
    for (; i < sizeP; i++) {
        resultat += P[i]*(pow(x, i));
    }
    return resultat;
}//nbCycles ~= 174368 donc ILP ~= 3/174368

float Ex9_Horner (float resultat, std::vector<float> P, int x, int i,  int sizeP) {
        for (; i > 0 ; i--) {
        if (i == sizeP) {
            resultat = P[i];
        }
        else {
            resultat = resultat*x + P[i];
        }
    }
    return resultat;
}//nbCycles ~= 47584 donc ILP ~= 3/47584

double & Ex10_reduce_mult_Qu1(double* V, size_t n) {
    __m256d x_vec, x_vec_V;
    double* temp = new double[4];

    x_vec_V = _mm256_set1_pd(1.0);

    //std::cout << V[0]*V[1]*V[2]*V[3]*V[4]*V[5] << std::endl;

    for (size_t i = 0; i < n-3; i+=4) {
        x_vec = _mm256_loadu_pd(V+i);
        x_vec_V = _mm256_mul_pd (x_vec_V, x_vec);
    }

    _mm256_storeu_pd (temp, x_vec_V);

    size_t mod = n -n %4;
    for (size_t j = n-1; j >= mod; j--) {
        temp[0] *= V[j];
    }

    for(size_t i = 1; i < 4; i++) {
        temp[0] *= temp[i];
    }

    return temp[0];
}

/*double & Ex10_reduce_mult_Qu2(int32_t V, size_t n) {
    __m256d x_vec, x_vec_V;
    int32_t* temp = new double[4];

    x_vec_V = _mm256_set1_pd(1.0);

    //std::cout << V[0]*V[1]*V[2]*V[3]*V[4]*V[5] << std::endl;

    for (size_t i = 0; i < n-3; i+=4) {
        x_vec = _mm256_loadu_pd(V+i);
        x_vec_V = _mm256_mul_pd (x_vec_V, x_vec);
    }

    _mm256_storeu_pd (temp, x_vec_V);

    size_t mod = n -n %4;
    for (size_t j = n-1; j >= mod; j--) {
        temp[0] *= V[j];
    }

    for(size_t i = 1; i < 4; i++) {
        temp[0] *= temp[i];
    }

    return temp[0];
}*/


double & Ex10_reduce_add_Qu1(double* V, size_t n) {
    __m256d x_vec, x_vec_V;
    double* temp = new double[4];

    x_vec_V = _mm256_set1_pd(0.0);

    std::cout << V[0]+V[1]+V[2]+V[3]+V[4]+V[5] << std::endl;

    for (size_t i = 0; i < n-3; i+=4) {
        x_vec = _mm256_loadu_pd(V+i);
        x_vec_V = _mm256_add_pd (x_vec_V, x_vec);
    }

    _mm256_storeu_pd (temp, x_vec_V);

    size_t mod = n -n %4;
    for (size_t j = n-1; j >= mod; j--) {
        temp[0] += V[j];
    }

    for(size_t i = 1; i < 4; i++) {
        temp[0] += temp[i];
    }

    return temp[0];

}

int32_t & Ex10_reduce_add_Qu2(int32_t* V, size_t n) {
    __m256i x_vec, x_vec_V;
    int32_t* temp = new int32_t[4];

    x_vec_V = _mm256_set1_epi32(0);

    std::cout << V[0]+V[1]+V[2]+V[3]+V[4]+V[5] << std::endl;

    for (size_t i = 0; i < n-3; i+=4) {
        x_vec = _mm256_loadu_si256((__m256 *)V+i);
        x_vec_V = _mm256_add_epi32(x_vec_V, x_vec);
    }

    _mm256_storeu_epi32(temp, x_vec_V);

    size_t mod = n -n %4;
    for (size_t j = n-1; j >= mod; j--) {
        temp[0] += V[j];
    }

    for(size_t i = 1; i < 4; i++) {
        temp[0] += temp[i];
    }

    return temp[0];

}

int32_t reduce3(int32_t *tab, size_t n) {
    int32_t sum = 1;
    for (int i = 0; i < n; i++) {
        sum *= tab[i];
    }
    return sum;
}

int main() {
  int n, N;
  EvalPerf PE;

  //Ex5 petit tableau
  std::vector<int> tab_Ex5_1(4);
  tab_Ex5_1[0] = 1;
  tab_Ex5_1[1] = 3;
  tab_Ex5_1[2] = 10;
  tab_Ex5_1[3] = 2;
  int tailleTab_Ex5_1 = tab_Ex5_1.size();
  std::vector<int> tabFinal_Ex5_1(tailleTab_Ex5_1);

  //Ex5 grand tableau
  std::vector<int> tab_Ex5_2(100);
  for (int i = 0; i < 100; i++) {
    tab_Ex5_2[i] = 100;
  }
  int tailleTab_Ex5_2 = tab_Ex5_2.size();
  std::vector<int> tabFinal_Ex5_2(tailleTab_Ex5_2);


  //Ex6 petit tableau int
  std::vector<float> tab_Ex6_1(4);
  tab_Ex6_1[0] = -5;
  tab_Ex6_1[1] = 3;
  tab_Ex6_1[2] = -7;
  tab_Ex6_1[3] = 4;
  int x_Ex6_1 = 2;
  float resultat_Ex6_1;

  //Ex6 grand tableau int
  std::vector<float> tab_Ex6_2(1000);
  for (int i = 0; i < 1000; i++) {
    tab_Ex6_2[i] = 30*i;
  }
  int x_Ex6_2 = 20;
  float resultat_Ex6_2;

  //Ex6 petit tableau float
  std::vector<float> tab_Ex6_3(10);
  tab_Ex6_3[0] = 10.6;
  tab_Ex6_3[1] = 5.2;
  tab_Ex6_3[2] = 34.7;
  tab_Ex6_3[3] = 2.1;
  tab_Ex6_3[4] = 7.9;
  tab_Ex6_3[6] = 16.9;
  tab_Ex6_3[7] = 65.6;
  tab_Ex6_3[8] = 12.8;
  tab_Ex6_3[9] = 7.4;
  int x_Ex6_3 = 3.5;
  float resultat_Ex6_3;

  //Ex6 grand tableau float
  std::vector<float> tab_Ex6_4(1000);
  for (int i = 0; i < 1000; i++) {
    tab_Ex6_4[i] = 45.6*i;
  }
  int x_Ex6_4 = 25.9;
  float resultat_Ex6_4;


  //Ex7
  std::vector<std::vector<float>> V = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  
  std::vector<float> resMult(3);
  resMult[0]=1;
  resMult[1]=1;
  resMult[2]=1;

  std::vector<float> resPlus(3);
  resPlus[0]=0;
  resPlus[1]=0;
  resPlus[2]=0;


  //Ex8
  std::vector<float> X(10000);
  std::vector<float> Y(10000);
  std::vector<float> Z(10000);

  for (int i = 0; i < 10000; i++) {
        X[i] = rand() % 100 + 1;
        Y[i] = rand() % 100 + 1;
        Z[i] = rand() % 100 + 1;
   }


  //Ex9 grand tableau int
  std::vector<float> tab_Ex9(1000);
  for (int i = 0; i < 1000; i++) {
    tab_Ex9[i] = 30*i;
  }
  int x_Ex9 = 20;
  float resultat_Ex9 = 0;
  int indice_Ex9_1 = 0;
  int sizeP_Ex9 = tab_Ex9.size();

  int indice_Ex9_2 = tab_Ex9.size()-1;


  //Ex10
  double* V_Ex10 = new double[6];
  V_Ex10[0] = 3;//1;
  V_Ex10[1] = 5;//2;
  V_Ex10[2] = 1;//3;
  V_Ex10[3] = 9;//4;
  V_Ex10[4] = 2;//5;
  V_Ex10[5] = 8;//6;
  double res_Ex10_mult = 1.;
  double res_Ex10_add = 0.;
  size_t n_Ex10 = 6;


  PE.start();
  
  /*Ex4_test();

  Ex5_sommePrefixe(tab_Ex5_1, tabFinal_Ex5_1);
  Ex5_sommePrefixe(tab_Ex5_2, tabFinal_Ex5_2);


  Ex6_puissances_alpha(resultat_Ex6_1, tab_Ex6_1, x_Ex6_1);
  Ex6_Horner(resultat_Ex6_1, tab_Ex6_1, x_Ex6_1);
  Ex6_puissances_alpha(resultat_Ex6_2, tab_Ex6_2, x_Ex6_2);
  Ex6_Horner(resultat_Ex6_2, tab_Ex6_2, x_Ex6_2);
  Ex6_puissances_alpha(resultat_Ex6_3, tab_Ex6_3, x_Ex6_3);
  Ex6_Horner(resultat_Ex6_3, tab_Ex6_3, x_Ex6_3);
  Ex6_puissances_alpha(resultat_Ex6_4, tab_Ex6_4, x_Ex6_4);
  Ex6_Horner(resultat_Ex6_4, tab_Ex6_4, x_Ex6_4);


  Ex7_reduce_mult(V, resMult);
  Ex7_reduce_plus(V, resPlus);


  //Ex8
  slowperformance1(X, Y, Z, 10000);
  slowperformance2(X, Y, Z, 10000);


  Ex9_puissances_alpha(resultat_Ex9, tab_Ex9, x_Ex9, indice_Ex9_1, sizeP_Ex9);
  Ex9_Horner (resultat_Ex9, tab_Ex9, x_Ex9, indice_Ex9_2, sizeP_Ex9);*/

  //double resultat1 = Ex10_reduce_mult_Qu1(V_Ex10, n_Ex10);
  double resultat2 = Ex10_reduce_add_Qu1(V_Ex10, n_Ex10);
  std::cout << resultat2 << std::endl;

  //pour récup un pointeur &(v[0]) ou v.data()

  PE.stop();

  std::cout << "nbr cycles: " << PE.cycles() << std::endl;
  std::cout << "nbr secondes: " << PE.seconds() << std::endl;
  std::cout << "nbr millisecondes: " << PE.millisecond() << std::endl;
  std::cout << "CPI = " << PE.CPI(N) << std::endl;
  std::cout << "IPC = " << PE.IPC(N) << std::endl;

  return 0;
}

