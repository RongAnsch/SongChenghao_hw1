#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    if (a.rows != b.rows || a.cols != b.cols)
    {
        printf("Error: Matrices dimensions do not match for addition.\n");
        return create_matrix(0, 0); // Return an empty matrix on error
    }

    Matrix result = create_matrix(a.rows, a.cols);
    int i=0;
    for (i; i < a.rows; ++i)
    {
        int j=0;
        for (j; j < a.cols; ++j)
        {
            result.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return result;
    return create_matrix(0, 0);
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    if (a.rows != b.rows || a.cols != b.cols)
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }

    Matrix result = create_matrix(a.rows, a.cols);
    int i=0;
    for (i; i < a.rows; ++i)
    {
        int j=0;
        for (j; j < a.cols; ++j)
        {
            result.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }
    return result;
    return create_matrix(0, 0);
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    if (a.cols != b.rows)
    {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }

    Matrix result = create_matrix(a.rows, b.cols);
    int i=0;
    for (i; i < a.rows; ++i)
    {
        int j=0;
        for (j; j < b.cols; ++j)
        {
            result.data[i][j] = 0.0;
            int k=0;
            for (k; k < a.cols; ++k)
            {
                result.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }
    return result;
    return create_matrix(0, 0);
}

Matrix scale_matrix(Matrix a, double k)
{
    Matrix result = create_matrix(a.rows, a.cols);
    int i=0;
    for (i; i < a.rows; ++i)
    {
        int j=0;
        for (j; j < a.cols; ++j)
        {
            result.data[i][j] = a.data[i][j] * k;
        }
    }
    return result;
    return create_matrix(0, 0);
}

Matrix transpose_matrix(Matrix a)
{
    Matrix result = create_matrix(a.cols, a.rows);
    int i=0;
    for (i; i < a.rows; ++i)
    {
        int j=0;
        for (j; j < a.cols; ++j)
        {
            result.data[j][i] = a.data[i][j];
        }
    }
    return result;
    return create_matrix(0, 0);
}

double det_matrix(Matrix a)
{
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0.0;
    }

    if (a.rows == 1)
    {
        return a.data[0][0];
    }

    if (a.rows == 2)
    {
        return a.data[0][0] * a.data[1][1] - a.data[0][1] * a.data[1][0];
    }

    double determinant = 0.0;
    int col=0;
    for (col; col < a.cols; ++col)
    {
        Matrix submatrix = create_matrix(a.rows - 1, a.cols - 1);
        int sub_i = 0, sub_j = 0;
        int i=1;
        for (i; i < a.rows; ++i)
        {
            sub_j = 0;
            int j=0;
            for (j; j < a.cols; ++j)
            {
                if (j == col)
                    continue;
                submatrix.data[sub_i][sub_j] = a.data[i][j];
                sub_j++;
            }
            sub_i++;
        }
        determinant += (col % 2 == 0 ? 1 : -1) * a.data[0][col] * det_matrix(submatrix);
    }
    return determinant;
}

Matrix inv_matrix(Matrix a)
{
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }

    double det = det_matrix(a);
    if (det == 0)
    {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }

    Matrix adj = create_matrix(a.rows, a.cols);
    int i=0;
    for (i; i < a.rows; ++i)
    {
        int j=0;
        for (j; j < a.cols; ++j)
        {
            Matrix submatrix = create_matrix(a.rows - 1, a.cols - 1);
            int sub_i = 0, sub_j = 0;
            int x=0;
            for (x; x < a.rows; ++x)
            {
                if (x == i)
                    continue;
                sub_j = 0;
                int y=0;
                for (y ; y < a.cols; ++y)
                {
                    if (y == j)
                        continue;
                    submatrix.data[sub_i][sub_j] = a.data[x][y];
                    sub_j++;
                }
                sub_i++;
            }
            adj.data[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * det_matrix(submatrix);
        }
    }

    Matrix inv = create_matrix(a.rows, a.cols);
    
    for (i = 0 ; i < a.rows; ++i)
    {
        int j=0;
        for (j ; j < a.cols; ++j)
        {
            inv.data[i][j] = adj.data[i][j] / det;
        }
    }

    return inv;
}

int rank_matrix(Matrix a)
{
    int rank = a.cols;
    Matrix temp = a;
    int row=0;
    for (row ; row < rank; row++)
    {
        if (temp.data[row][row] != 0)
        {
            int col=0;
            for (col ; col < a.rows; col++)
            {
                if (col != row)
                {
                    double mult = temp.data[col][row] / temp.data[row][row];
                    int i=0;
                    for (i ; i < rank; i++)
                    {
                        temp.data[col][i] -= mult * temp.data[row][i];
                    }
                }
            }
        }
        else
        {
            int reduce = 1;
            int i=row+1;
            for (i; i < a.rows; i++)
            {
                if (temp.data[i][row] != 0)
                {
                    int j=0;
                    for (j ; j < rank; j++)
                    {
                        double tmp = temp.data[row][j];
                        temp.data[row][j] = temp.data[i][j];
                        temp.data[i][j] = tmp;
                    }
                    reduce = 0;
                    break;
                }
            }
            if (reduce)
            {
                rank--;
                int i=0;
                for (i ; i < a.rows; i++)
                {
                    temp.data[i][row] = temp.data[i][rank];
                }
            }
            row--;
        }
    }
    return rank;
    return 0;
}

double trace_matrix(Matrix a)
{
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0.0;
    }

    double trace = 0.0;
    int i=0;
    for (i ; i < a.rows; ++i)
    {
        trace += a.data[i][i];
    }
    return trace;
    return 0;
}

void print_matrix(Matrix a)
{
    int i=0;
    for (i ; i < a.rows; i++)
    {
        int j=0;
        for (j ; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}