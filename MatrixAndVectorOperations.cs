using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PyroElectricsFEM
{
    class Matrix3
    {
        private double a11 = 0;
        private double a12 = 0;
        private double a13 = 0;
        private double a21 = 0;
        private double a22 = 0;
        private double a23 = 0;
        private double a31 = 0;
        private double a32 = 0;
        private double a33 = 0;

        public Matrix3() { }
        public Matrix3(double a, double b, double c, double a1, double b1, double c1, double a2, double b2, double c2)
        {
            a11 = a;
            a12 = b;
            a13 = c;

            a21 = a1;
            a22 = b1;
            a23 = c1;

            a31 = a2;
            a32 = b2;
            a33 = c2;
        }
        public double A11
        {
            get { return a11; }
            set { a11 = value; }
        }
        public double A12
        {
            get { return a12; }
            set { a12 = value; }
        }
        public double A13
        {
            get { return a13; }
            set { a13 = value; }
        }
        public double A21
        {
            get { return a21; }
            set { a21 = value; }
        }
        public double A22
        {
            get { return a22; }
            set { a22 = value; }
        }
        public double A23
        {
            get { return a23; }
            set { a23 = value; }
        }
        public double A31
        {
            get { return a31; }
            set { a31 = value; }
        }
        public double A32
        {
            get { return a32; }
            set { a32 = value; }
        }
        public double A33
        {
            get { return a33; }
            set { a33 = value; }
        }
        public Matrix3 InvertedMatrix()
        {
            double det = a11*a22*a33-a11*a23*a32-a12*a21*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31;
            if (det == 0)
                return null;
            else
            {
                double k11 = (a22*a33-a23*a32) / det;
                double k12 = -(a21*a33-a23*a31) / det;
                double k13 = (a21*a32-a22*a31) / det;
                double k21 = -(a12*a33-a13*a32) / det;
                double k22 = (a11*a33-a13*a31) / det;
                double k23 = -(a11*a32-a12*a31) / det;
                double k31 = (a12*a23-a13*a22) / det;
                double k32 = -(a11*a23-a13*a21) / det;
                double k33 = (a11*a22-a12*a21) / det;

                return new Matrix3(k11,k21,k31,k12,k22,k32,k13,k23,k33);
            }
        }
        public static Matrix3 operator *(Matrix3 A, Matrix3 B)
        {
            Matrix3 newMatrix = new Matrix3();

            newMatrix.A11 = A.A11 * B.A11 + A.A12 * B.A21 + A.A13 * B.A31;
            newMatrix.A12 = A.A11 * B.A12 + A.A12 * B.A22 + A.A13 * B.A32;
            newMatrix.A13 = A.A11 * B.A13 + A.A12 * B.A23 + A.A13 * B.A33;

            newMatrix.A21 = A.A21 * B.A11 + A.A22 * B.A21 + A.A23 * B.A31;
            newMatrix.A22 = A.A21 * B.A12 + A.A22 * B.A22 + A.A23 * B.A32;
            newMatrix.A23 = A.A21 * B.A13 + A.A22 * B.A23 + A.A23 * B.A33;

            newMatrix.A31 = A.A31 * B.A11 + A.A32 * B.A21 + A.A33 * B.A31;
            newMatrix.A32 = A.A31 * B.A12 + A.A32 * B.A22 + A.A33 * B.A32;
            newMatrix.A33 = A.A31 * B.A13 + A.A32 * B.A23 + A.A33 * B.A33;

            return newMatrix;
        }
        public static Vector3 operator *(Matrix3 A, Vector3 B)
        {
            Vector3 newVector = new Vector3();
            newVector.A1 = A.A11 * B.A1 + A.A12 * B.A2+A.A13*B.A3;
            newVector.A2 = A.A21 * B.A1 + A.A22 * B.A2+A.A23*B.A3;
            newVector.A3 = A.A31 * B.A1 + A.A32 * B.A2+A.A33*B.A3;
            return newVector;
        }
        public static Matrix3 operator +(Matrix3 A, Matrix3 B)
        {
            Matrix3 newMatrix = new Matrix3();

            newMatrix.A11 = A.A11 + B.A11;
            newMatrix.A12 = A.A12 + B.A12;
            newMatrix.A13 = A.A13 + B.A13;

            newMatrix.A21 = A.A21 + B.A21;
            newMatrix.A22 = A.A22 + B.A22;
            newMatrix.A23 = A.A23 + B.A23;

            newMatrix.A31 = A.A31 + B.A31;
            newMatrix.A32 = A.A32 + B.A32;
            newMatrix.A33 = A.A33 + B.A33;

            return newMatrix;
        }
        public static Matrix3 operator -(Matrix3 A, Matrix3 B)
        {
            Matrix3 newMatrix = new Matrix3();

            newMatrix.A11 = A.A11 - B.A11;
            newMatrix.A12 = A.A12 - B.A12;
            newMatrix.A13 = A.A13 - B.A13;

            newMatrix.A21 = A.A21 - B.A21;
            newMatrix.A22 = A.A22 - B.A22;
            newMatrix.A23 = A.A23 - B.A23;

            newMatrix.A31 = A.A31 - B.A31;
            newMatrix.A32 = A.A32 - B.A32;
            newMatrix.A33 = A.A33 - B.A33;

            return newMatrix;
        }
        public static Matrix3 operator -(Matrix3 A)
        {
            Matrix3 newMatrix = new Matrix3();

            newMatrix.A11 = -A.A11;
            newMatrix.A12 = -A.A12;
            newMatrix.A13 = -A.A13;

            newMatrix.A21 = -A.A21;
            newMatrix.A22 = -A.A22;
            newMatrix.A23 = -A.A23;

            newMatrix.A31 = -A.A31;
            newMatrix.A32 = -A.A32;
            newMatrix.A33 = -A.A33;

            return newMatrix;
        }
        public override string ToString()
        {
            return $"{a11},{a12},{a13}\n{a21},{a22},{a23}\n{a31},{a32},{a33}\n";
        }
        public object Clone()
        {
            return new Matrix3(this.a11, this.a12, this.a13, this.a21, this.a22, this.a23, this.a31, this.a32, this.a33);
        }

    }
    class Vector3
    {
        private double a1 = 0;
        private double a2 = 0;
        private double a3 = 0;

        public Vector3() { }
        public Vector3(double a, double b, double c)
        {
            a1 = a;
            a2 = b;
            a3 = c;
        }
        public double A1
        {
            get { return a1; }
            set { a1 = value; }
        }
        public double A2
        {
            get { return a2; }
            set { a2 = value; }
        }
        public double A3
        {
            get { return a3; }
            set { a3 = value; }
        }
        public static Vector3 operator *(Vector3 B, Matrix3 A)
        {
            Vector3 newVector = new Vector3();

            newVector.A1 = A.A11 * B.A1 + A.A12 * B.A2+A.A13*B.A3;
            newVector.A2 = A.A21 * B.A1 + A.A22 * B.A2+A.A23*B.A3;
            ///????
            newVector.A3= A.A31 * B.A1 + A.A32 * B.A2 + A.A33 * B.A3;

            return newVector;
        }
        public static Vector3 operator +(Vector3 A, Vector3 B)
        {
            Vector3 newVector = new Vector3();

            newVector.A1 = A.A1 + B.A1;
            newVector.A2 = A.A2 + B.A2;
            newVector.A3 = A.A3 + B.A3;

            return newVector;
        }
        public static Vector3 operator -(Vector3 A, Vector3 B)
        {
            Vector3 newVector = new Vector3();

            newVector.A1 = A.A1 - B.A1;
            newVector.A2 = A.A2 - B.A2;
            newVector.A3 = A.A3 - B.A3;

            return newVector;
        }
        public static Vector3 operator -(Vector3 A)
        {
            Vector3 newVector = new Vector3();

            newVector.A1 = -A.A1;
            newVector.A2 = -A.A2;
            newVector.A3 = -A.A3;
            return newVector;
        }
        public override string ToString()
        {
            return $"{a1}, {a2}, {a3}";
        }
        public object Clone()
        {
            return new Vector3(this.a1, this.a2, this.a3);
        }
    }
}
