#include <bits/stdc++.h>
using namespace std;
const int N = 3e6 + 10;
const int M = 2e7 + 10;

int n, m, Num_n; 
int row_ptr[M], col_idx[M]; 
double val[M];
double b[N];
double x_n[N]; 

void prepare_read(ifstream &file1, ifstream &file2)
{
    int h, k;
    file1 >> n >> m >> h;  
    for (int i = 0; i <= n; i++) file1 >> row_ptr[i];
    file1 >> Num_n;
    for (int i = 0; i < Num_n; i++) file1 >> col_idx[i];
    file1 >> k;
    for (int i = 0; i < Num_n; i++) 
    {
        file1 >> val[i];
    }
    for (int i = 0; i < n; i++) 
    {
        string str;
        file2 >> str;
        b[i] = stod(str);
    }
}

void sparse_matrix_vector_multiply(double *vec, double *result)
{
    for (int i = 0; i < n; i++)
    {
        result[i] = 0;
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++)
        {
            result[i] += val[j] * vec[col_idx[j]];
        }
    }
}


double dot_product(double *vec1, double *vec2)
{
    long double result = 0; 
    for (int i = 0; i < n; i++)
    {
        result += (long double)vec1[i] * vec2[i];
    }
    return (double)result;
}

void vector_add(double *vec1, double *vec2, double *result, double alpha)
{
    for (int i = 0; i < n; i++)
    {
        result[i] = vec1[i] + alpha * vec2[i];
    }
}

void vector_subtract(double *vec1, double *vec2, double *result, double alpha)
{
    for (int i = 0; i < n; i++)
    {
        result[i] = vec1[i] - alpha * vec2[i];
    }
}
double r[N], p[N], Ap[N];

void solve_conjugate_gradient()
{
    
    for (int i = 0; i < n; i++) x_n[i] = 0; 

    sparse_matrix_vector_multiply(x_n, r);
    vector_subtract(b, r, r, 1.0); // r = b - Ax
    for (int i = 0; i < n; i++) p[i] = r[i]; // p = r

    double rs_old = dot_product(r, r);
    int Iteration = 0;
    double Tolerance = 1e-7;

    while (sqrt(rs_old) > Tolerance && Iteration < 2000)
    {
        sparse_matrix_vector_multiply(p, Ap); // Ap = A * p
        double dot_pAp = dot_product(p, Ap);

        if (fabs(dot_pAp) < 1e-12) {
            for (int i = 0; i < n; i++) {
                val[row_ptr[i]] += 1e-8; 
            }
            continue; 
        }
        double alpha = rs_old / dot_pAp;

        vector_add(x_n, p, x_n, alpha); // x = x + alpha * p
        vector_subtract(r, Ap, r, alpha); // r = r - alpha * Ap

        double rs_new = dot_product(r, r);

        if (sqrt(rs_new) <= Tolerance) break;

        double beta = rs_new / rs_old;
        vector_add(r, p, p, beta); // p = r + beta * p

        rs_old = rs_new;
        Iteration++;

        cout << "Iterations: " << Iteration << "    Error:  " << sqrt(rs_new) << endl;
    }
        cout << "Final Solution Vector (x):" << endl;
    for (int i = 0; i < n; i++) {
        cout << x_n[i] << " ";
    }
    cout << endl;
}

int main()
{
    ifstream file1("1.in");
    ifstream file2("2.in");
    freopen("Conjugate Gradient Method.out","w",stdout);

    prepare_read(file1, file2); 
    solve_conjugate_gradient();

    file1.close();
    file2.close();

    return 0;
}
