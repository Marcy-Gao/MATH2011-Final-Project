#include <bits/stdc++.h>
using namespace std;
const int N = 3e6 + 10;
const int M = 2e7 + 10;

int n, m, Num_n; // n * m matrix; Num_n: nonzero number
int row_ptr[M], col_idx[M]; 
double val[M], x_n[N];
double b[N];
double res1[N]; // value of D^(-1)*b
double dig[N]; //diagonal element
double Tolerance;
double res2[N];
inline void read(int &x)  // quick input to reduce time 
{
    x = 0; int w = 1;
    char ch = getchar();
    while (!isdigit(ch)) {if (ch == '-') w = -1; ch = getchar();}
    while (isdigit(ch)) {x = x * 10 + (ch ^ 48); ch = getchar();}
    x = x * w;
}

void prepare_read(ifstream &file1, ifstream &file2)
{
    int h, k;
    file1 >> n >> m >> h;  
    for (int i = 0; i <= n; i++) file1 >> row_ptr[i];
    file1 >> Num_n;
    for (int i = 0; i < Num_n; i++) file1 >> col_idx[i];
    file1 >> k;
    for (int i = 0; i < Num_n; i++) file1 >> val[i];
    for (int i = 0; i < n; i++) 
    {
        string str;
        file2 >> str;
        b[i] = stod(str);
        // cout << b[i] << endl;
    }
    cout << "!!" << endl;
}


void work() // Calculate the next x;
{
    
    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = row_ptr[i] + 1; j < row_ptr[i + 1]; j++)
        {
            if (col_idx[j] == i) continue;
            sum -= val[j] * x_n[col_idx[j]];
        }
        x_n[i] = (sum + b[i]) / val[row_ptr[i]] ; //use the new value instead of last x_n
    }
}

double error_calculate()
{
    for (int i = 0; i < n; i++)
    {
        res2[i] = 0;
        res2[i] += val[row_ptr[i]] * x_n[i];
        for (int j = row_ptr[i] + 1; j < row_ptr[i + 1]; j++)
        {
            res2[i] += val[j] * x_n[col_idx[j]];
        }

    }  
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        res2[i] -= b[i];
        double t = res2[i];
        sum += t * t;
    } 
    double error = sqrt(sum);
    return error;
}

void output( int Iteration, double error)
{
    cout << "Final Answer:  " << error << endl;
    for (int i = 0; i < n; i++) cout << x_n[i] << " ";
    cout << endl;
}

void solve_iteration()
{

    int Iteration = 0;
    Tolerance = 1e-8; double error;
    for (int i = 0; i < n; i++) x_n[i] = 0;
    while (1)
    {
        // last[N] = x_n[N];
        work();
        error = error_calculate();
        Iteration++;
        if (error <= Tolerance || Iteration > 2000) 
        {
            break;
        }
        
        cout << "Iterations: " << Iteration << "    Error:  " << error << endl;
    }
    output( Iteration, error);
}

int main()
{
    ifstream file1("1.in");
    ifstream file2("2.in");
    freopen("Gauss-Seidel Method.out","w",stdout);
    if (!file1.is_open() || !file2.is_open()) {
        cout << "Error opening files!" << endl;
        return -1;
    }

    prepare_read(file1, file2); 
 
    solve_iteration();

    file1.close();
    file2.close();

    return 0;
}
