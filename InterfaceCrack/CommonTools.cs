using System;

namespace InterfaceCrack
{
    public static class CommonTools
    {
        public static double Integrate(Func<double, double> func, double b)
        {
            double sum = 0;
            double h = b / 100;

            double x = h;

            while (x < b)
            {
                sum += 4 * func(x);
                x += h;
                if (x >= b) break;
                sum += 2 * func(x);
                x += h;
            }
            sum = (h / 3) * (sum + func(0) - func(b));

            return sum;
        }

        public static void CalculateInfinity(Func<double, double>[] functions, out double[] _infinityQ, out double[] _infinityT)
        {
            _infinityQ = new double[4];
            _infinityT = new double[4];
            for (int i = 0; i < functions.Length; i++)
            {
                int step = 0;
                double Sn = 0, Sn1 = 100;
                while (Math.Abs(Sn1 - Sn) > 0.000001)
                {
                    step++;
                    Sn = Sn1;
                    Sn1 = functions[i](step);
                }
                _infinityQ[i] = Sn1;
                _infinityT[i] = step;
            }
        }
    }
}
