#include <stdio.h>
#include <math.h>

double v(double r) {
  delta = 2.8; // a.u.
  r_c = 6.7; //a.u.
  if(r<=(r_c-0.5*delta)) && (r>=(r_c+0.5*delta){
    return -0.30
  }  
  else{
   return 0; 
  }
}

double f(double r, double E, int l) {
    double rsqr = r * r;
    double M = 1, Z = 1;
    double a1 = 16.039, a2 = 2.007, a3 = -25.543, a4 = 4.525, a5 = 0.961, a6 = 0.443;
    double dum2 = (Z + a1 * exp(-a2 * r) + a3 * r * exp(-a4 * r) + a5 * exp(-a6 * r)) / r;
    return (v(r) - E) * 2 * M + (l * (l + 1) / rsqr) - (2 * M * dum2);
}

int main() {
    int i, l, n, preci = 0;
    double h = 0.003, d = 0.1, r[60000], u[60000], w[60000], R1[2000], R2[2000];
    double U[60000], R[60000], r1, r2, U1, U2, E = -0.8, pie = 3.14;
    double check = 0, dum, dum2;
    FILE *fp = fopen("final_bwfn_r_14au.in", "w");
    FILE *fp2 = fopen("unormbwfn_total.dat", "w");

    if (!fp || !fp2) {
        printf("Error opening file!\n");
        return 1;
    }

    l = 1;
    for (n = 0; n <= 1000; n++) {
        E += d;
        r[0] = 0.001;
        r[1] = r[0] + h;
        u[0] = pow(r[0], l + 1);
        u[1] = pow(r[1], l + 1);

        w[0] = u[0] * (1 - (h * h * f(r[0], E, l)) / 12);
        w[1] = u[1] * (1 - (h * h * f(r[1], E, l)) / 12);

        for (i = 2; i <= 8340; i++) {
            r[i] = r[i - 1] + h;
            w[i] = -w[i - 2] + (2.0 * w[i - 1]) + (h * h * f(r[i - 1], E, l) * u[i - 1]);
            u[i] = w[i] / (1 - (h * h * f(r[i], E, l) / 12.0));
        }

        for (i = 0; i <= 8340; i++) {
            U[i] = u[i];
            R[i] = U[i] / r[i];
        }

        r1 = r[4704];
        r2 = r[4740];
        U1 = U[4704];
        U2 = U[4740];
        R1[n] = R[4704];
        R2[n] = R[4740];

        if (n >= 1) {
            check = R1[n] * R1[n - 1];
            if (check < 0) {
                preci++;
                E -= d;
                d /= 10;
                n = 0;
            }
        }

        if (preci > 8) {
            printf("Energy is %lf\n", E);
            printf("Precision is %d\n", preci);
            break;
        }

        fprintf(fp, "%lf  %d  %lf   %lf\n", E, n, R1[n], check);
    }

    for (i = 0; i <= 8340; i++) {
        fprintf(fp2, "%lf  %lf   %lf  %lf\n", E, r[i], U[i], R[i]);
    }

    fclose(fp);
    fclose(fp2);
    return 0;
}

