using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using static PyroElectricsFEM.MathHelper;

namespace PyroElectricsFEM
{
    public partial class Form1 : Form
    {
        #region fields
        double L;
        int N;
        double[] xi;
        double step;
        Function[] fi;
        Function[] fiDerivatives;
        double ro;
        double sigma;
        double alpha;
        double w=0;
        double f=0;
        double c;
        double E;
        double D;
        double g;
        double hShtrih;
        double lambda;
        double pSmall;
        double[,] CMatrix;
        double[,] EMatrix;
        double[,] YMatrix;
        double[,] GMatrix;
        double[,] KMatrix;
        double[,] PMatrix;
        double[] LV;
        double[] RQ;
        double[] UXI;
        #endregion
        public Form1()
        {
            InitializeComponent();
        }
        public double[,] Fill1Koef(double koeficient)
        {
            double[,] CMatrix = new double[N + 1, 3];

            CMatrix[0, 0] = 0;
            CMatrix[0, 1] = IntegrateGauss((x) => (koeficient * fiDerivatives[0](x) * fiDerivatives[0](x)), xi[0], xi[1]);
            CMatrix[0, 2] = IntegrateGauss((x) => (koeficient * fiDerivatives[0](x) * fiDerivatives[1](x)), xi[0], xi[1]);

            for (int i = 1; i < N; i++)
            {
                CMatrix[i, 0] = IntegrateGauss((x) => (koeficient * fiDerivatives[i - 1](x) * fiDerivatives[i](x)), xi[i - 1], xi[i]);
                CMatrix[i, 1] = IntegrateGauss((x) => (koeficient * fiDerivatives[i](x) * fiDerivatives[i](x)), xi[i - 1], xi[i]);
                CMatrix[i, 1] += IntegrateGauss((x) => (koeficient * fiDerivatives[i](x) * fiDerivatives[i](x)), xi[i], xi[i + 1]);
                CMatrix[i, 2] = IntegrateGauss((x) => (koeficient * fiDerivatives[i](x) * fiDerivatives[i + 1](x)), xi[i], xi[i + 1]);

            }

            CMatrix[N, 0] = IntegrateGauss((x) => (koeficient * fiDerivatives[N - 1](x) * fiDerivatives[N](x)), xi[N - 1], xi[N]);
            CMatrix[N, 1] = IntegrateGauss((x) => (koeficient * fiDerivatives[N](x) * fiDerivatives[N](x)), xi[N - 1], xi[N]);
            CMatrix[N, 2] = 0;
            double[,] CMatrixNew = new double[N + 1, N + 1];
            for (int i = 0; i <= N; i++)
            {
                CMatrixNew[i, i] = CMatrix[i, 1];
                if (i > 0)
                {
                    CMatrixNew[i, i - 1] = CMatrix[i, 0];
                }
                if (i < N)
                {
                    CMatrixNew[i, i + 1] = CMatrix[i, 2];
                }
            }
            return CMatrixNew;
        }
        public double[,] FillYMatrix(double c, double alpha) 
        {
            double[,] CMatrix = new double[N + 1, 3];

            CMatrix[0, 0] = 0;
            CMatrix[0, 1] = IntegrateGauss((x) => (c*alpha * fi[0](x) * fiDerivatives[0](x)), xi[0], xi[1]);
            CMatrix[0, 2] = IntegrateGauss((x) => (c * alpha * fi[0](x) * fiDerivatives[1](x)), xi[0], xi[1]);

            for (int i = 1; i < N; i++)
            {
                CMatrix[i, 0] = IntegrateGauss((x) => (c * alpha * fi[i - 1](x) * fiDerivatives[i](x)), xi[i - 1], xi[i]);
                CMatrix[i, 1] = IntegrateGauss((x) => (c * alpha * fi[i](x) * fiDerivatives[i](x)), xi[i - 1], xi[i]);
                CMatrix[i, 1] += IntegrateGauss((x) => (c * alpha * fi[i](x) * fiDerivatives[i](x)), xi[i], xi[i + 1]);
                CMatrix[i, 2] = IntegrateGauss((x) => (c * alpha * fi[i](x) * fiDerivatives[i + 1](x)), xi[i], xi[i + 1]);

            }

            CMatrix[N, 0] = IntegrateGauss((x) => (c * alpha * fi[N - 1](x) * fiDerivatives[N](x)), xi[N - 1], xi[N]);
            CMatrix[N, 1] = IntegrateGauss((x) => (c * alpha * fi[N](x) * fiDerivatives[N](x)), xi[N - 1], xi[N]);
            CMatrix[N, 2] = 0;
            double[,] CMatrixNew = new double[N + 1, N + 1];
            for (int i = 0; i <= N; i++)
            {
                CMatrixNew[i, i] = CMatrix[i, 1];
                if (i > 0)
                {
                    CMatrixNew[i, i - 1] = CMatrix[i, 0];
                }
                if (i < N)
                {
                    CMatrixNew[i, i + 1] = CMatrix[i, 2];
                }
            }
            return CMatrixNew;
        }
        public double[,] FillPMatrix(double p)
        {
            double[,] CMatrix = new double[N + 1, 3];

            CMatrix[0, 0] = 0;
            CMatrix[0, 1] = IntegrateGauss((x) => (p * fi[0](x) * fiDerivatives[0](x)), xi[0], xi[1]);
            CMatrix[0, 2] = IntegrateGauss((x) => (p * fi[0](x) * fiDerivatives[1](x)), xi[0], xi[1]);

            for (int i = 1; i < N; i++)
            {
                CMatrix[i, 0] = IntegrateGauss((x) => (p * fi[i - 1](x) * fiDerivatives[i](x)), xi[i - 1], xi[i]);
                CMatrix[i, 1] = IntegrateGauss((x) => (p * fi[i](x) * fiDerivatives[i](x)), xi[i - 1], xi[i]);
                CMatrix[i, 1] += IntegrateGauss((x) => (p * fi[i](x) * fiDerivatives[i](x)), xi[i], xi[i + 1]);
                CMatrix[i, 2] = IntegrateGauss((x) => (p * fi[i](x) * fiDerivatives[i + 1](x)), xi[i], xi[i + 1]);

            }

            CMatrix[N, 0] = IntegrateGauss((x) => (p * fi[N - 1](x) * fiDerivatives[N](x)), xi[N - 1], xi[N]);
            CMatrix[N, 1] = IntegrateGauss((x) => (p * fi[N](x) * fiDerivatives[N](x)), xi[N - 1], xi[N]);
            CMatrix[N, 2] = 0;
            double[,] CMatrixNew = new double[N + 1, N + 1];
            for (int i = 0; i <= N; i++)
            {
                CMatrixNew[i, i] = CMatrix[i, 1];
                if (i > 0)
                {
                    CMatrixNew[i, i - 1] = CMatrix[i, 0];
                }
                if (i < N)
                {
                    CMatrixNew[i, i + 1] = CMatrix[i, 2];
                }
            }
            return CMatrixNew;
        }

        private void TestBtn_Click(object sender, EventArgs e)
        {
            vChart.Series[0].Points.Clear();
            qChart.Series[0].Points.Clear();
            zChart.Series[0].Points.Clear();

            L = double.Parse(LTextBox.Text);
            N = Convert.ToInt32(NTextBox.Text);
            c = Convert.ToDouble(cXTextBox.Text);
            E = Convert.ToDouble(eXTextBox.Text);
            alpha = Convert.ToDouble(alphaXTextBox.Text);
            g = Convert.ToDouble(gXTextBox.Text);
            pSmall = Convert.ToDouble(pXTextBox.Text);
            lambda = Convert.ToDouble(lambdaXTextBox.Text);
            ro = Convert.ToDouble(roXTextBox.Text);
            sigma = Convert.ToDouble(sigmaTextBox.Text);
            D = Convert.ToDouble(DTextBox.Text);
            hShtrih = Convert.ToDouble(hLTextBox.Text);


            step = L / N;
            //courant functions
            fi = new Function[N + 1];
            fiDerivatives = new Function[N + 1];
            for (int i = 0; i <=N; i++)
            {
                CourantFunctions c = new CourantFunctions(L,i,N);
                fi[i] = new Function(c.fi);
                fiDerivatives[i] = new Function(c.fiDerivative);
            }

            //xi forming
            xi = new double[N + 1];
            xi[0] = 0;
            for (int i = 1; i <= N; i++)
            {
                xi[i] = xi[i - 1] + step;
            }

            foreach (var x in xi)
            {
                Console.WriteLine(x);
            }
            //forming Cmatrix
            CMatrix = Fill1Koef(c);
            for (int i = 0; i <=N; i++)
            {
                for (int j = 0; j <=N; j++)
                {
                    Console.Write($"{CMatrix[i, j]}   ");
                }
                Console.WriteLine();
            }
            //forming Ematrix
            EMatrix = Fill1Koef(E);
            //Ymatrix
            YMatrix = FillYMatrix(c,alpha);
            for (int i = 0; i <= N; i++)
            {
                for (int j = 0; j <= N; j++)
                {
                    Console.Write($"{YMatrix[i,j]}  ");
                }
                Console.WriteLine(  );
            }
            //GMatrix
            GMatrix = Fill1Koef(g);
            //KMatrix
            KMatrix = Fill1Koef(lambda);
            //PMatrix
            PMatrix = Fill1Koef(pSmall);

            double[,] MainMatrix = new double[3 * (N + 1), 3 * (N + 1)];

            for (int i = 0; i <=N; i++)
            {
                for (int j = 0; j <= N; j++)
                {
                    //first row
                    MainMatrix[i, j] = CMatrix[i, j];
                    MainMatrix[i, j + N + 1] = EMatrix[i, j];
                    MainMatrix[i, j + 2*(N + 1)] = (-1)*YMatrix[i, j];
                    //second row
                    MainMatrix[i + N + 1, j] = EMatrix[i, j];
                    MainMatrix[i + N + 1, j + N + 1] = (-1) * GMatrix[i, j];
                    MainMatrix[i + N + 1, j + 2 * (N + 1)] = PMatrix[i, j];
                    //third row 
                    MainMatrix[i + 2 * (N + 1), j + 2 * (N + 1)] = KMatrix[i, j];

                }
            }
            MainMatrix[0, 0] = Math.Pow(10, 20);
            MainMatrix[N+1, N+1] = Math.Pow(10, 20);
            MainMatrix[2*N + 2, 2 * N + 2] = Math.Pow(10, 20);

            for (int i = 0; i < MainMatrix.GetLength(0); i++)
            {
                for (int j = 0; j < MainMatrix.GetLength(0); j++)
                {
                    Console.Write($"{MainMatrix[i, j]}   ");
                }
                Console.WriteLine();
            }

            Console.WriteLine("----------------------------------------------------------");
            LV = new double[N + 1];
            for (int i = 0; i <= N; i++)
            {
                LV[i] = IntegrateGauss((x) => (ro * f * fi[i](x)), 0, L) + sigma * fi[i](L);
            }

            RQ = new double[N + 1];
            for (int i = 0; i <= N; i++)
            {
                RQ[i] = D*fi[i](L);
            }

            UXI = new double[N + 1];
            for (int i = 0; i <= N; i++)
            {
                UXI[i] = IntegrateGauss((x) => (ro * w * fi[i](x)), 0, L) - hShtrih * fi[i](L);
            }

            for (int i = 0; i <= N; i++)
            {
                Console.WriteLine($"LV {i} = {LV[i]}");
                Console.WriteLine($"RQ {i} = {RQ[i]}");
                Console.WriteLine($"UXI {i} = {UXI[i]}");
            }

            var MainVectorTmp= LV.Concat(RQ).Concat(UXI).ToArray();
            double[] MainVector = new double[3 * (N + 1)];
            for (int i = 0; i < MainVectorTmp.Length; i++)
            {
                MainVector[i] = MainVectorTmp[i];
            }
            double[] result = Gauss(MainMatrix, MainVector);
            foreach (var item in result)
            {
                Console.WriteLine(item);
            }
            double[] uKoefs = new double[N + 1];
            double[] pKoefs = new double[N + 1];
            double[] ThettaKoefs = new double[N + 1];
            for (int i = 0; i <=N ; i++)
            {
                uKoefs[i] = result[i];
                pKoefs[i] = result[N + 1 + i];
                ThettaKoefs[i] = result[2 * N + 2 + i];
            }
            for (int i = 0; i <=N; i++)
            {
                vChart.Series[0].Points.AddXY(xi[i], uKoefs[i]);
                qChart.Series[0].Points.AddXY(xi[i], pKoefs[i]);
                zChart.Series[0].Points.AddXY(xi[i], ThettaKoefs[i]);
            }

            //double[,] testM = { {1,1,1 }, {2,1,3 }, {7,8,9 } };
            //double[] testV = { 3,2,1};
            //double[] testResult = Gauss(testM, testV);
            //foreach (var item in testResult)
            //{
            //    Console.WriteLine(item);
            //}

        }
    }
}
