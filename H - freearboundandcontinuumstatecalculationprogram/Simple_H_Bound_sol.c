#include <stdio.h>
#include <math.h>

double v(double r);
double f(double r, double E, int l);

int main()
{
    int i, l, n, preci = 0;
    double h = 0.003, d = 0.1, r[60000], u[60000], w[60000];
    double U[60000], R[60000], r1, r2, U1, U2, R1[2000], R2[2000], E, check = 0;
    double fr, rsqr, Z = 1, pie = 3.14, M = 1;
    FILE *fp, *fp2;

    fp = fopen("final_bwfn_r_14au.in", "w");
    fp2 = fopen("unormbwfn_total.dat", "w");

    if (!fp || !fp2) {
        printf("Error opening files!\n");
        return 1;
    }

    // 1s Electron Binding Energy of Hydrogen
    E = -0.8;  // Corrected from -0.7
    l = 0;      // Hydrogen ground state

    for (n = 0; n <= 1000; n++) {
        E += d;  // Increment energy step
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
            printf("Final Energy: %lf \n", E);
            printf("Precision Level: %d\n", preci);
            break;
        }

        fprintf(fp, "%lf  %d  %lf   %lf \n", E, n, R1[n], check);
    }

    for (i = 0; i <= 8340; i++) {
        fprintf(fp2, "%lf  %lf   %lf  %lf \n", E, r[i], U[i], R[i]);
    }

    fclose(fp);
    fclose(fp2);
    return 0;
}

// Coulomb potential function for Hydrogen
double v(double r) {
    return -1.0 / r;  // Hydrogen potential (-Z/r)
}

// Effective potential function
double f(double r, double E, int l) {
    double rsqr = r * r, M = 1, Z = 1;
    return (v(r) - E) * 2 * M + (l * (l + 1) / rsqr) - (2 * M * v(r));
}

