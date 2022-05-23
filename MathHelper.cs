using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PyroElectricsFEM
{
    class MathHelper
    {
        
        public delegate double Function(double x);
        public static double IntegrateGauss(Function f, double a, double b)
        {
            double[] x = new double[8];
            x[0] = -0.96028986;
            x[7] = -x[0];
            x[1] = -0.79666648;
            x[6] = -x[1];
            x[2] = -0.52553242;
            x[5] = -x[2];
            x[3] = -0.18343464;
            x[4] = -x[3];

            double[] A = new double[8];
            A[0] = A[7] = 0.10122854;
            A[1] = A[6] = 0.22238103;
            A[2] = A[5] = 0.31370664;
            A[3] = A[4] = 0.36268378;

            double[] t = new double[8];
            for (int i = 0; i < 8; i++)
            {
                t[i] = (a + b) / 2 + ((b - a) / 2) * x[i];
            }

            double integral = 0;
            for (int i = 0; i < 8; i++)
            {
                integral += (b - a) * (A[i] * f(t[i])) / 2;

            }

            return integral;
        }
        public static double Derivative(Function f, double x)
        {
            double eps = Math.Pow(10, 20);
            return ((f(x + eps) - f(x)) / eps);

        }
        public static List<double> ThomasMethod(List<double> mainDiagonal, List<double> upperDiagonal, List<double> lowerDiagonal, List<double> bVector)
        {
            int dimensionSize = mainDiagonal.Count;
            List<double> cShtrih = Enumerable.Repeat(0.0, dimensionSize).ToList();
            List<double> dShtrih = Enumerable.Repeat(0.0, dimensionSize).ToList();
            List<double> resultVector = Enumerable.Repeat(0.0, dimensionSize).ToList();

            //koeficients calculation
            cShtrih[0] = upperDiagonal[0] / mainDiagonal[0];
            dShtrih[0] = Convert.ToDouble(bVector[0] / mainDiagonal[0]);


            for (int i = 1; i < dimensionSize - 1; i++)
            {
                cShtrih[i] = upperDiagonal[i] / (mainDiagonal[i] - lowerDiagonal[i] * cShtrih[i - 1]);
                dShtrih[i] = (bVector[i] - lowerDiagonal[i] * dShtrih[i - 1]) / (mainDiagonal[i] - lowerDiagonal[i] * cShtrih[i - 1]);
            }
            dShtrih[dimensionSize - 1] = (bVector[dimensionSize - 1] - lowerDiagonal[dimensionSize - 1] * dShtrih[dimensionSize - 2]) / (mainDiagonal[dimensionSize - 1] - lowerDiagonal[dimensionSize - 1] * cShtrih[dimensionSize - 2]);

            //result calculation
            resultVector[dimensionSize - 1] = dShtrih[dimensionSize - 1];
            for (int i = dimensionSize - 2; i >= 0; i--)
            {
                resultVector[i] = dShtrih[i] - cShtrih[i] * resultVector[i + 1];
            }

            return resultVector;
        }
        public static double[] Multiply(double[,] a, double[] v)
        {
            double[] res = new double[v.Length];

            if (a.GetLength(0) == v.Length && a.GetLength(1) == 3)
            {
                res[0] = a[0, 1] * v[0] + a[0, 2] * v[1];
                for (int i = 1; i < v.Length - 1; i++)
                {
                    res[i] = a[i, 0] * v[i - 1] + a[i, 1] * v[i] + a[i, 2] * v[i + 1];

                }

                res[v.Length - 1] = a[v.Length - 1, 0] * v[v.Length - 2] + a[v.Length - 1, 1] * v[v.Length - 1];

            }

            return res;
        }
        public static double[] Gauss(double[,] A, double[] B)
        {
            int N = B.Length;
            int i, j, k;
            double R;

            //Прямий хід
            for (i = 0; i <= N - 2; i++)
            {
                k = i; // номер рядка, який будемо міняти з і-тим

                R = Math.Abs(A[i, i]);
                for (j = i + 1; j <= N - 1; j++)
                    if (Math.Abs(A[j, i]) >= R)
                    {
                        k = j;
                        R = Math.Abs(A[j, i]);
                    }


                // якщо треба, то міняємо рядки матриці та відповідні елементи в стовпчику
                if (k != i)
                {
                    R = B[k];
                    B[k] = B[i];
                    B[i] = R;

                    for (j = i; j < N; j++)
                    {
                        R = A[k, j];
                        A[k, j] = A[i, j];
                        A[i, j] = R;
                    }
                }


                R = A[i, i];
                B[i] = B[i] / R;
                for (j = 1; j <= N - 1; j++)
                    A[i, j] = A[i, j] / R;

                for (k = i + 1; k <= N - 1; k++)
                {
                    R = A[k, i];
                    B[k] = B[k] - R * B[i];

                    A[k, i] = 0;
                    for (j = i + 1; j <= N - 1; j++)
                        A[k, j] = A[k, j] - R * A[i, j];
                }
            }



            //Обернений хід

            double[] X = new double[N];
            X[N - 1] = B[N - 1] / A[N - 1, N - 1];
            for (i = N - 2; i >= 0; i--)
            {
                R = B[i];
                for (j = i + 1; j <= N - 1; j++)
                    R -= A[i, j] * X[j];

                X[i] = R;
            }

            return X;
        }
        // Метод прогонки розв'язування системи Ax = F, де А - блочно-тридіагональна матриця
        public static Vector3[] TridiagonalMatrixMethod(Matrix3[,] A, Vector3[] F)
        {
            int N = F.Length;

            //коефіцієнти прогонки
            Matrix3[] a = new Matrix3[N];
            Vector3[] b = new Vector3[N];

            a[1] = -A[0, 2] * A[0, 1].InvertedMatrix();
            b[1] = F[0] * A[0, 1].InvertedMatrix();

            for (int i = 2; i < N; ++i)
            {
                a[i] = -A[i - 1, 2] * ((A[i - 1, 0] * a[i - 1] + A[i - 1, 1]).InvertedMatrix());
                b[i] = (F[i - 1] - A[i - 1, 0] * b[i - 1]) * ((A[i - 1, 0] * a[i - 1] + A[i - 1, 1]).InvertedMatrix());
            }

            Vector3[] X = new Vector3[N];

            X[N - 1] = (F[N - 1] - A[N - 1, 0] * b[N - 1]) * (A[N - 1, 1] + A[N - 1, 0] * a[N - 1]).InvertedMatrix();
            for (int i = N - 2; i >= 0; --i)
                X[i] = a[i + 1] * X[i + 1] + b[i + 1];

            return X;

        }

        //public static Vector3[] BlockTridiagonalGauss(Matrix3[,] A, Vector3[] F)
        //{

        //    int N = F.Length;

        //    Vector3[] res = new Vector3[N];

        //    //Прямий хід
        //    Matrix3 restMatrix = new Matrix3();
        //    Vector3 restVector = new Vector3();
        //    for (int i = 0; i < N - 1; i++)
        //    {
        //        A[i, 1] = A[i, 1] + restMatrix;
        //        F[i] = F[i] + restVector;

        //        Matrix3 A11 = (Matrix2)A[i, 1].Clone();
        //        Matrix2 A12 = (Matrix2)A[i, 2].Clone();
        //        Matrix2 A21 = (Matrix2)A[i + 1, 0].Clone();

        //        Vector2 V1 = (Vector2)F[i].Clone();
        //        Vector2 V2 = (Vector2)F[i + 1].Clone();

        //        A[i, 1] = new Matrix2(1, 0, 0, 1);
        //        A[i + 1, 0] = new Matrix2();
        //        A[i, 2] = A11.InvertedMatrix() * A12;
        //        restMatrix = (-A21 * (A11.InvertedMatrix())) * A12;

        //        F[i] = A11.InvertedMatrix() * V1;
        //        restVector = (-A21 * (A11.InvertedMatrix())) * V1;


        //    }
        //    A[N - 1, 1] = A[N - 1, 1] + restMatrix;
        //    F[N - 1] = F[N - 1] + restVector;

        //    //Обернений хід
        //    res[N - 1] = (A[N - 1, 1].InvertedMatrix()) * F[N - 1];
        //    for (int i = N - 2; i >= 0; --i)
        //    {
        //        res[i] = F[i] - A[i, 2] * res[i + 1];

        //    }

        //    return res;



        //}
    }
}
