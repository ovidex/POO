#include <iostream>

class Complex{
private:
    double real;
    double imaginary;
public:

    Complex(){
        real = 0;
        imaginary = 0;
    }

    Complex(int a, int b){
        real = a;
        imaginary = b;
    }

    ~Complex(){

    }

    friend Complex& operator+(Complex, Complex);
    friend Complex& operator*(Complex, Complex);
    friend Complex& operator/(Complex, Complex);
    friend std::ostream& operator<<(std::ostream&, Complex);
    friend std::istream& operator>>(std::istream&, Complex&);

    friend class ComplexMatrix;
};

Complex& operator+(Complex a, Complex b){
    Complex *c = new Complex;
    c->real = a.real + b.real;
    c->imaginary = a.imaginary + b.imaginary;

    return *c;
}

Complex& operator*(Complex a, Complex b){
    Complex *c = new Complex;
    double re = a.real;
    double im = a.imaginary;

    c->real = (re * b.real) - (im * b.imaginary);
    c->imaginary = (re * b.imaginary) + (im * b.real);

    return *c;
}

Complex& operator/(Complex a, Complex b){
    Complex *c = new Complex;

    double div = b.real * b.real + b.imaginary * b.imaginary;
    b.imaginary = - b.imaginary;

    Complex upp = a * b;

    c->real = upp.real/div;
    c->imaginary = upp.imaginary/div;

    return *c;
}

std::ostream& operator<<(std::ostream& out, Complex a){
    out<<a.real;
    if(a.imaginary >= 0)
        out<<"+"<<a.imaginary<<"i";
    else
        out<<a.imaginary<<"i";

    return out;
}

std::istream& operator>>(std::istream& in, Complex& a){
    char sign, im;
    in>>a.real>>a.imaginary>>im;

    return in;
}

class ComplexMatrix{
private:
    Complex bp[100][100];
    int rows, cols;

    void cofactor(Complex mat[100][100], Complex tmp[100][100], int p, int q, int n){
        int i = 0, j =0 ;

        for(int x = 0; x < rows; x++){
            for(int y = 0; y < cols; y++){
                if(x != p && y != q){
                    tmp[i][j++] = mat[x][y];

                    if(j == rows - 1){
                        j = 0;
                        i++;
                    }
                }
            }
        }
    }

    Complex determinant(Complex mat[100][100], int n){
        Complex result;

        if(n == 1)
            return mat[0][0];

        Complex tmp[100][100];

        Complex sign(1, 0);

        for(int f = 0; f < n; f++){
            cofactor(mat, tmp, 0, f, n);
            Complex aux = mat[0][f] * determinant(tmp, n - 1);
            aux = sign * aux;
            result = result + aux;

            sign.real = -sign.real;
        }

        return result;
    }

    void adjoint(Complex A[100][100], Complex adj[100][100]){
        Complex sign(1, 0), temp[100][100];

        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                cofactor(A, temp, i, j, rows);

                sign.real = ((i+j)%2==0)? 1: -1;

                adj[j][i] = (sign)*(determinant(temp, rows-1));
            }
        }
    }

public:
    ComplexMatrix(int n, int m){
        rows = n;
        cols = m;
    }

    ~ComplexMatrix(){

    }

    Complex getDeterminant(){
        if(rows != cols){
            exit(1);
        }
        return determinant(bp, rows);
    }

    ComplexMatrix& inverse(){

        Complex det = getDeterminant();

        if(det.real == 0)
            exit(1);

        Complex adj[100][100];

        adjoint(this->bp, adj);

        ComplexMatrix* ret = new ComplexMatrix(rows, cols);

        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                ret->bp[i][j] = adj[i][j]/det;

        return *ret;

    }

    friend ComplexMatrix& operator+(ComplexMatrix, ComplexMatrix);
    friend ComplexMatrix& operator*(ComplexMatrix, ComplexMatrix);
    friend std::istream& operator>>(std::istream&, ComplexMatrix&);
    friend std::ostream& operator<<(std::ostream&, ComplexMatrix&);
};

ComplexMatrix& operator+(ComplexMatrix a, ComplexMatrix b){
    if(a.rows != b.rows || a.cols != b.cols){
        exit(1);
    }

    ComplexMatrix *c = new ComplexMatrix(a.rows, a.cols);

    for(int i = 0; i < a.rows; i++)
        for(int j = 0; j < a.cols; j++)
            c->bp[i][j] = a.bp[i][j] + b.bp[i][j];

    return *c;
}

ComplexMatrix& operator*(ComplexMatrix a, ComplexMatrix b){
    if(a.rows != b.cols || a.cols != b.rows)
        exit(1);

    ComplexMatrix *c = new ComplexMatrix(a.rows, b.cols);

    for(int i = 0; i < a.rows; i++){
        for(int j = 0; j < b.cols; j++){
            for(int k = 0; k < a.cols; k++){
                c->bp[i][j] = c->bp[i][j] + (a.bp[i][k] * b.bp[k][j]);
            }
        }
    }

    return *c;
}

std::istream& operator>>(std::istream& in, ComplexMatrix& m){
    for(int i = 0; i < m.rows; i++)
        for(int j = 0; j < m.cols; j++)
            in>>m.bp[i][j];

    return in;
}

std::ostream& operator<<(std::ostream& out, ComplexMatrix &m){
    for(int i = 0; i < m.rows; i++){
        for(int j = 0; j < m.cols; j++)
            out<<m.bp[i][j]<<"\t";
        out<<"\n";
    }

    return out;
}

int main(){
    ComplexMatrix a(2, 2);
    std::cin>>a;
    ComplexMatrix b(2, 2);
    b = a.inverse();
    std::cout<<b;
    return 0;
}
