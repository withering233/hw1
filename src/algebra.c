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
    Matrix c = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
            c.data[i][j] = a.data[i][j] + b.data[i][j];
    }
    return c;
    if (a.rows != b.rows || a.cols != b.cols)
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    Matrix c = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
            c.data[i][j] = a.data[i][j] - b.data[i][j];
    }
    return c;
    if (a.rows != b.rows || a.cols != b.cols)
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    Matrix c = create_matrix(a.rows, b.cols);
    // printf("%d %d",a.cols,b.rows);
    if (a.cols != b.rows)
    {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < b.cols; j++)
        {
            c.data[i][j] = 0;
            for (int k = 0; k < a.cols; k++)
            {
                c.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }
    return c;
}

Matrix scale_matrix(Matrix a, double k)
{
    Matrix c = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
            c.data[i][j] = a.data[i][j] * k;
    }
    return c;
}

Matrix transpose_matrix(Matrix a)
{
    Matrix c = create_matrix(a.cols, a.rows);
    for (int i = 0; i < a.cols; i++)
    {
        for (int j = 0; j < a.rows; j++)
            c.data[i][j] = a.data[j][i];
    }
    return c;
}

double det(double matrix[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int order);
double minor(double matrix[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int r, int c, int n);

double det_matrix(Matrix a)
{
    double result;
    if (a.rows != a.cols)
    {
        printf("Error: Matrix a must be a square matrix.\n");
        return 0;
    }

    result = det(a.data, a.rows);

    return result;
}

double det(double matrix[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int order)
{
    double result = 0.0, sign = 1, cofactor;
    int i;

    if (order == 1)
    {
        result = matrix[0][0];
    }
    else
    {
        for (i = 0; i < order; i++)
        {
            cofactor = minor(matrix, i, 0, order);
            result += sign * matrix[i][0] * cofactor;
            sign *= -1;
        }
    }

    return result;
}

double minor(double matrix[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int r, int c, int n)
{
    int original_i, original_j, i, j;
    double result = 0.0;
    double matrix2[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            original_i = i;
            original_j = j;
            if (i == r || j == c)
                ;
            else
            {
                if (i > r)
                    i--;
                if (j > c)
                    j--;
                matrix2[i][j] = matrix[original_i][original_j];
                i = original_i;
                j = original_j;
            }
        }
    if (n >= 2)
        result = det(matrix2, n - 1);
    return result;
}

Matrix inv_matrix(Matrix a)
{
    Matrix b = create_matrix(a.rows, a.cols);
    if (a.rows != a.cols || det_matrix(a) == 0)
    {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }

    double augmentedMatrix[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE * 2];
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.rows; j++)
        {
            augmentedMatrix[i][j] = a.data[i][j];
            augmentedMatrix[i][j + a.rows] = (i == j) ? 1.0 : 0.0;
        }
    }
    for (int i = 0; i < a.rows; i++)
    {
        double divisor = augmentedMatrix[i][i];
        for (int j = 0; j < 2 * a.rows; j++)
        {
            augmentedMatrix[i][j] /= divisor;
        }
        for (int k = 0; k < a.rows; k++)
        {
            if (k != i)
            {
                double factor = augmentedMatrix[k][i];
                for (int j = 0; j < 2 * a.rows; j++)
                {
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                }
            }
        }
    }
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.rows; j++)
        {
            b.data[i][j] = augmentedMatrix[i][j + a.rows];
        }
    }
    return b;
}

int rank_matrix(Matrix a)
{
    int rank; 
    int i, j, k, lead;
    if (a.rows <= a.cols)
        rank = a.rows;
    else
        rank = a.cols;

    for (i = 0; i < a.rows; i++)
    {
        lead = 0;
        while (lead < a.cols && a.data[i][lead] == 0)
        {
            lead++;
        }
        if (lead == a.cols)
        {
            rank--;
        }
        else
        {
            for (j = i + 1; j < a.rows; j++)
            {
                float scalar = a.data[j][lead] / a.data[i][lead];
                for (k = lead; k < a.cols; k++)
                {
                    a.data[j][k] -= scalar * a.data[i][k];
                }
            }
        }
    }

    return rank;
}

double trace_matrix(Matrix a)
{
    if (a.rows != a.cols)
    {
        printf("Error: Matrix a must be a square matrix.\n");
        return 0;
    }
    double sum = 0.0;
    for (int i = 0; i < a.rows; i++)
    {
        sum += a.data[i][i];
    }
    return sum;
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}