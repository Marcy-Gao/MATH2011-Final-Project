#include <bits/stdc++.h>
#include <sys/time.h>
using namespace std;

const int N = 3e6 + 10;  // Maximum number of rows
const int M = 2e7 + 10;  // Maximum number of non-zero elements in matrix

int n, m, Num_n;         // n * m matrix; Num_n: non-zero elements count
int row_ptr[M], col_idx[M]; 
double val[M];           // Values of the non-zero elements
double b[N];             // Right-hand side vector
double res1[N], res2[N], res3[N];  // intermediate results
double x_n[N];           // Current solution vector
double dig[N];           // Diagonal elements (not used currently)
long double Tolerance;   // Convergence tolerance

inline void read(int &x) {
    x = 0; int w = 1;
    char ch = getchar();
    while (!isdigit(ch)) {if (ch == '-') w = -1; ch = getchar();}
    while (isdigit(ch)) {x = x * 10 + (ch ^ 48); ch = getchar();}
    x = x * w;
}

void prepare_read(ifstream &file1, ifstream &file2) {
    int h, k;
    file1 >> n >> m >> h;  
    for (int i = 0; i <= n; i++) file1 >> row_ptr[i];
    file1 >> Num_n;
    for (int i = 0; i < Num_n; i++) file1 >> col_idx[i];
    file1 >> k;
    for (int i = 0; i < Num_n; i++) file1 >> val[i];
    for (int i = 0; i < n; i++) {
        string str;
        file2 >> str;
        b[i] = stod(str);
    }
}

void work() {
    // Calculate the next approximation of x_n
    for (int i = 0; i < n; i++) {
        res3[i] = 0;
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            if (col_idx[j] != i) {  // Skip diagonal elements
                res3[i] += val[j] * res2[col_idx[j]];
            }
        }
        res3[i] += val[row_ptr[i]] * res2[i]; // Adding the diagonal term
    }

    long double sum1 = 0, sum2 = 0; 
    for (int i = 0; i < n; i++) {
        sum1 += res2[i] * res2[i];
        sum2 += res2[i] * res3[i];
    }

    long double alpha = sum1 / sum2;
    for (int i = 0; i < n; i++) {
        x_n[i] -= alpha * res2[i];
    }
}

long double error_calculate() {
    // Calculate the error ||Ax - b||
    // fill(res2, res2 + n, 0); // Reset res2
    for (int i = 0; i < n; i++) 
    {
        res2[i] = 0;
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            if (col_idx[j] != i) {  // Skip diagonal elements
                res2[i] += val[j] * x_n[col_idx[j]];
            }
        }
    
        res2[i] += val[row_ptr[i]] * x_n[i]; // Adding the diagonal term
    }
    for (int i = 0; i < n; i++) {
        // cout << i << " " << res2[i] << endl;
        res2[i] += b[i]; // res2 = Ax - b
    }

    long double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += res2[i] * res2[i];
    }
    return sqrt(sum);  // Return the norm ||Ax - b||
}

void output( int Iteration, double error)
{
    cout << "Final Answer:  " << error << endl;
    for (int i = 0; i < n; i++) cout << x_n[i] << " ";
    cout << endl;
}

void solve_iteration() {
    int Iteration = 0;
    Tolerance = 1e-8; // Convergence tolerance
    long double error = 0;

    fill(x_n, x_n + n, 0);  // Initial guess for x_n = 1
    copy(b, b + n, res2);      // res2 = b initially
    while (1) {
        work();  // Update x_n using the steepest descent method
        error = error_calculate();  // Calculate error ||Ax - b||
        Iteration++;
        
        if (error <= Tolerance || Iteration > 2000) {
            break;
        }
        cout << "Iterations: " << Iteration << "    Error:  " << error << endl;
    }
    output(Iteration, error);
}

int main() {
    ifstream file1("1.in");
    ifstream file2("2.in");
    freopen("Steepest Descent Method.out","w",stdout);


    prepare_read(file1, file2); 
    solve_iteration();

    file1.close();
    file2.close();

    return 0;
}
