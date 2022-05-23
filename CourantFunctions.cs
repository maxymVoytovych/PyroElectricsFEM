using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PyroElectricsFEM
{
    class CourantFunctions
    {
        private int i;
        private int NumberOfFiniteElements;
        private double Length;

        public CourantFunctions() { i = 0; NumberOfFiniteElements = 10; Length = 1.0; }
        public CourantFunctions(int number, int total) { i = number; NumberOfFiniteElements = total; Length = 1.0; }
        public CourantFunctions(double length, int functionId, int numberOfElements) { Length = length; i = functionId; NumberOfFiniteElements = numberOfElements; }
        public double fi(double x)
        {
            double h = Length / NumberOfFiniteElements;

            double f;

            if (x < 0 || x > Length)
                f = 0;
            else
            {

                if ((i - 1) * h < x && x <= i * h)
                {
                    f = (x - (i - 1) * h) / h;
                }
                else
                    if (i * h < x && x <= (i + 1) * h)
                {
                    f = ((i + 1) * h - x) / h;
                }
                else
                    f = 0;
            }

            return f;
        }
        public double fiDerivative(double x)
        {
            double h = Length / NumberOfFiniteElements;

            double f;

            if (x < 0 || x > Length)
                f = 0;
            else
            {

                if ((i - 1) * h < x && x <= i * h)
                {
                    f = NumberOfFiniteElements / Length;
                }
                else
                    if (i * h < x && x <= (i + 1) * h)
                {
                    f = -NumberOfFiniteElements / Length;
                }
                else
                    f = 0;
            }

            return f;
        }
    }
}
