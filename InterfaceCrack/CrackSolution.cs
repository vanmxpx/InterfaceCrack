using System;
using System.ComponentModel;
using System.Globalization;
using System.IO;
using System.Linq;

namespace InterfaceCrack
{
    public class CrackSolution
    {
        private NumberFormatInfo format;
        private Action<string> printer;
        private BackgroundWorker backgroundWorker;
        private Matrix _matrixA;

        const int N = 100;
        double p11, p22;
        double nu1, nu2;
        double mu1, mu2, mu, b1, b2;
        double s, p, h, a, g;
        double alpha11, alpha12, alpha21, alpha22;
        double alpha1b, alpha2b, alpha1c, alpha1d, alpha2c, alpha2d;
        double delta1, delta2, a11, a12, a13, a21, a22, a23, a31, a32, a33, a41, a42, a43;
        Func<double, double> r11, r21, r22, r12, r31, r32, r41, r42, r51, r52, r61, r62;
        Func<double, double> n11, n12, n21, n22;
        Func<double, double> q11, q12, q21, q22;
        Func<double, double, double> K11, K12, K21, K22;
        Func<double, double> TETTA, deltaMain;
        double[] _infinityT = new double[4];
        double[] _infinityQ = new double[4];
        double[] _sArray;
        double[] _tArray;
        double[] p11array;
        double[] detarray;
        double[] Answer;

        public CrackSolution(Action<string> printer, BackgroundWorker backgroundWorker)
        {
            this.backgroundWorker = backgroundWorker;
            this.printer = printer;
            p11array = new double[N];
            detarray = new double[N];

            InvokeFunctions();
        }

        void InvokeFunctions()
        {
            format = new NumberFormatInfo();
            format.NumberGroupSeparator = ",";
            format.NumberDecimalSeparator = ".";

            r11 = t => (alpha1c * s + 2 * alpha1d * alpha11 * alpha11) * Math.Exp(-2 * t * h * alpha11) / delta1;
            r12 = t => 2 * alpha1d * s * Math.Exp(-t * h * (alpha11 + alpha12)) / delta1;
            r21 = t => -4 * alpha1c * alpha11 * alpha11 * Math.Exp(-t * h * (alpha11 + alpha12)) / delta1;
            r22 = t => -(alpha1c * s + 2 * alpha1d * alpha11 * alpha11) * Math.Exp(-2 * t * h * alpha12) / delta1;
            r31 = t => (a11 * r11(t) + a12 * r21(t) + a13) / delta2;
            r32 = t => (a21 * r12(t) + a22 * r22(t) + a23) / delta2;
            r41 = t => (a31 * r11(t) + a32 * r21(t) + a33) / delta2;
            r42 = t => (a41 * r12(t) + a42 * r22(t) + a43) / delta2;
            r51 = t => a11 * r11(t) + a12 * r21(t) - a21 * r31(t) - a22 * r41(t) - a11;
            r52 = t => a11 * r12(t) + a12 * r22(t) - a21 * r32(t) - a22 * r42(t) - a12;
            r61 = t => a11 * a11 * r11(t) + r21(t) - a21 * a21 * r31(t) - r41(t) + a11 * a11;
            r62 = t => a11 * a11 * r12(t) + r22(t) - a21 * a21 * r32(t) - r42(t) + 1;

            n11 = t => alpha1c * r11(t) + alpha1d * r21(t) - alpha1c;
            n12 = t => alpha1c * r12(t) + alpha1d * r22(t) - alpha1d;
            n21 = t => 2 * alpha11 * alpha11 * r11(t) + s * r21(t) + 2 * alpha11 * alpha11;
            n22 = t => 2 * alpha11 * alpha11 * r12(t) + s * r22(t) + s;

            deltaMain = t => r51(t) * r62(t) - r61(t) * r52(t);

            q11 = t => (n11(t) * r62(t) - n12(t) * r61(t)) / deltaMain(t);
            q21 = t => (n21(t) * r62(t) - n22(t) * r61(t)) / deltaMain(t);
            q12 = t => (n12(t) * r51(t) - n11(t) * r52(t)) / deltaMain(t);
            q22 = t => (n22(t) * r51(t) - n21(t) * r52(t)) / deltaMain(t);

            K11 = (j, x) => CommonTools.Integrate(t => (q11(t) - _infinityQ[0]) * Math.Cos(t * j) * Math.Cos(t * x), _infinityT[0]) / _infinityQ[1];
            K12 = (j, x) => CommonTools.Integrate(t => (q12(t) - _infinityQ[1]) * Math.Sin(t * j) * Math.Cos(t * x), _infinityT[1]) / _infinityQ[1];
            K21 = (j, x) => CommonTools.Integrate(t => (q21(t) - _infinityQ[2]) * Math.Cos(t * j) * Math.Sin(t * x), _infinityT[2]) / -_infinityQ[2];
            K22 = (j, x) => CommonTools.Integrate(t => (q22(t) - _infinityQ[3]) * Math.Sin(t * j) * Math.Sin(t * x), _infinityT[3]) / -_infinityQ[2];

            TETTA = Si => - g * Math.Log(1 - Si * Si);
        }

        internal double[] CalculateWithNeuralNetworkTest(double h, double a, double v1, double v2, double mu1, double mu2)
        {
            var number = new[] { h, v1, v2 };
            double[][] w1 = new double[3][];
            double[][] w2 = new double[500][]; 
            double[][] wout = new double[500][];
            string[] w1lines = File.ReadAllLines(@"crack_wghts_layer1_4000ep.txt").ToArray();
            string[] w2lines = File.ReadAllLines(@"crack_wghts_layer2_4000ep.txt").ToArray();
            string[] woutlines = File.ReadAllLines(@"crack_wghts_out_4000ep.txt").ToArray();

            // разобрать в массив
            for (int i = 0; i < 3; i++)
            {
                double[] row = w1lines[i].Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Select(v => Double.Parse(v, format)).ToArray();
                w1[i] = row;
            }

            for (int i = 0; i < 500; i++)
            {
                double[] row = w2lines[i].Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Select(v => Double.Parse(v, format)).ToArray();
                w2[i] = row;
            }

            for (int i = 0; i < 500; i++)
            {
                double[] row = woutlines[i].Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Select(v => Double.Parse(v, format)).ToArray();
                wout[i] = row;
            }

            var inputs = number;
            double tmp;
            double[] outputs1 = new double[500];
            double[] outputs2 = new double[500];
            double[] result = new double[27];

            for (int row = 0; row < 500; row++)
            {
                //inputs sinapses
                tmp = 0;
                for (int col = 0; col < 3; col++)//(3 dendrids in 2000 neurons)
                    tmp += inputs[col] * w1[col][row];
                outputs1[row] = normalizeAnswer(tmp);
            }
            for (int row = 0; row < 500; row++)
            {
                //inputs sinapses
                tmp = 0;
                for (int col = 0; col < 500; col++)////(2000 dendrids in 2000 neurons)
                    tmp += outputs1[col] * w2[col][ row];
                outputs2[row] = normalizeAnswer(tmp);
            }
            //hidden sinapses
            for (int row = 0; row < 27; row++)
            {
                tmp = 0;
                for (int col = 0; col < 500; col++) //(2000 dendrids in 27 neurons)
                    tmp += outputs2[col] * wout[col][row];
                result[row] = normalizeAnswer(tmp);
            }

            printer("*Critical P: " + result[0] + Environment.NewLine);

            double axys = result[2] - (result[2] - result[1]);
            double[] res = new double[50];
            for (int i = 0; i < 26; i++)
            {
                if (i < 24)
                    res[26 - i + 23] = result[26 - i];
                res[i] = result[26 - i];
            }
            return res;
        }

        internal double[] CalculateWithNeuralNetwork(double h, double a, double v1, double v2, double mu1, double mu2)
        {
            return analize(new []{h, v1, v2});
        }

        double[] analize(double[] number)
        {
            double[][] w = new double[3][]; ;
            double[][] wout = new double[2000][];
            string[] wlines = File.ReadAllLines(@"crack_wghts_test.txt").ToArray();
            string[] woutlines = File.ReadAllLines(@"crack_wghts_out_test.txt").ToArray();

            // разобрать в массив
            for (int i = 0; i < 3; i++)
            {
                double[] row = wlines[i].Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Select(v => Double.Parse(v, format)).ToArray();
                w[i] = row;
            }

            for (int i = 0; i < 2000; i++)
            {
                double[] row = woutlines[i].Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Select(v => Double.Parse(v, format)).ToArray();
                wout[i] = row;
            }

            double[] result  = new double[27];

            var inputs = number;
            double tmp;
            double[] outputs = new double[2000];
            for (int row = 0; row < 2000; row++)
            {
                //inputs sinapses
                tmp = 0;
                for (int col = 0; col < 3; col++)
                    tmp += inputs[col] * w[col][row];
                outputs[row] = normalizeAnswer(tmp);
            }
            //hidden sinapses
            for (int row = 0; row < 27; row++)
            {
                tmp = 0;
                for (int col = 0; col < 2000; col++)
                    tmp += outputs[col] * wout[col][row];
                result[row] = normalizeAnswer(tmp);
            }
            printer("*Critical P: " + result[0] + Environment.NewLine);

            double[] res = new double[50];
            for (int i = 0; i < 26; i++)
            {
                if(i < 24)
                    res[26 - i + 23] = result[26 - i];
                res[i] = result[26 - i];
            }
            return res;
        }

        double normalizeAnswer(double answer)
        {
            return Math.Log(answer + Math.Sqrt(answer * answer + 1));
        }
        double[] globalSolutione = new double[46];
        private double div = 0.01, dec = 0;
        public double[] Calculate(int n, double h, double a, double nu1, double nu2, double mu1, double mu2) //коеф пуассона, модуль сдвига
        {
            //backgroundWorker.ReportProgress(0);

            this.h = h;
            this.a = a;

            _matrixA = new Matrix(2*n);
            _sArray = new double[n];
            _tArray = new double[n - 1];

            mu = mu2 / mu1;

            for (int i = 1; i < N; i++)
            {
                p11 = i * div + dec;

                p22 = (1 - nu1) * mu * p11 / (1 - nu2);

                alpha11 = Math.Sqrt(1 - p11 * (1 - 2 * nu1) / (2 * (1 - nu1)));
                alpha21 = Math.Sqrt(1 - p22 * (1 - 2 * nu2) / (2 * mu * (1 - nu2)));
                alpha12 = Math.Sqrt(1 - p11);
                alpha22 = Math.Sqrt(1 - p22 / mu);

                //alpha1a = 2 * ((1 - nu1) * alpha11 - nu1 * Math.Pow(alpha11, 3)) / (1 - 2 * nu1);
                //alpha2a = 2 * mu * ((1 - nu2) * alpha21 - nu2 * Math.Pow(alpha21, 3)) / (1 - 2 * nu2);
                alpha1b = 2 * alpha12;
                alpha2b = 2 * mu * alpha22;
                //2
                alpha1c = 2 * (nu1 * alpha11 - (1 - nu1) * Math.Pow(alpha11, 3)) / (1 - 2 * nu1);
                alpha2c = 2 * mu * (nu2 * alpha21 - (1 - nu2) * Math.Pow(alpha21, 3)) / (1 - 2 * nu2);
                alpha1d = -alpha1b;
                alpha2d = -alpha2b;

                s = 1 + alpha12 * alpha12; p = 1 + alpha22 * alpha22;

                delta1 = alpha1c * s - 2 * alpha1d * alpha11 * alpha11;
                delta2 = alpha2c * p - 2 * alpha2d * alpha21 * alpha21;

                a11 = alpha1c * p - 2 * alpha2d * alpha11 * alpha11 / mu;
                a12 = alpha1d * p - alpha2d * s / mu;
                a13 = -(alpha1c * p + 2 * alpha2d * alpha11 * alpha11 / mu);
                a21 = a11;
                a22 = a12;
                a23 = -(alpha1d * p + alpha2d * s / mu);
                //3
                a31 = 2 * alpha2c * alpha11 * alpha11 / mu - 2 * alpha1c * alpha21 * alpha21;
                a32 = alpha2c * s / mu - 2 * alpha1d * alpha21 * alpha21;
                a33 = 2 * alpha2c * alpha11 * alpha11 / mu + 2 * alpha1c * alpha21 * alpha21;
                a41 = a31;
                a42 = a32;
                a43 = alpha2c * s / mu + 2 * alpha1d * alpha21 * alpha21;
                var aii = new[] {a11, a12, a13, a21, a22, a23 , a31, a32, a33 , a41, a42, a43 };
                CommonTools.CalculateInfinity(new[] { q11, q12, q21, q22 }, out _infinityQ, out _infinityT);

                b1 = _infinityQ[0] / _infinityQ[1];
                b2 = _infinityQ[3] / _infinityQ[2];

                g = Math.Log((1 + Math.Sqrt(b1 * b2)) / (1 - Math.Sqrt(b1 * b2))) / (2 * Math.PI);

                #region Filling the matrix
                //Main functions
                double[] CosTetta = new double[n];
                double[] SinTetta = new double[n];
                _matrixA.Null();
                for (int k = 0; k < n; k++)
                {
                    _sArray[k] = Math.Cos((2 * k + 1) * Math.PI / (2 * n));
                    CosTetta[k] = Math.Cos(TETTA(_sArray[k]));
                    SinTetta[k] = Math.Sin(TETTA(_sArray[k]));
                    if (k < n - 1)
                        _tArray[k] = Math.Cos((k + 1) * Math.PI / n);
                }

                for (int k = 0; k < n - 1; k++)
                {   //Main diaghonal (1,n-1)
                    _matrixA.Body[k][k] += 0.5 * b1 * CosTetta[k] / Math.Sqrt(1 - Math.Pow(_sArray[k], 2));
                    _matrixA.Body[k][k + 1] += 0.5 * b1 * CosTetta[k + 1] / Math.Sqrt(1 - Math.Pow(_sArray[k + 1], 2));
                    _matrixA.Body[k][n + k] -= 0.5 * b1 * SinTetta[k] / Math.Sqrt(1 - Math.Pow(_sArray[k], 2));
                    _matrixA.Body[k][n + k + 1] -= 0.5 * b1 * SinTetta[k + 1] / Math.Sqrt(1 - Math.Pow(_sArray[k + 1], 2));
                    //Main diaghonal (n, 2n-2)
                    _matrixA.Body[n - 1 + k][k] -= 0.5 * b2 * SinTetta[k] / Math.Sqrt(1 - Math.Pow(_sArray[k], 2));
                    _matrixA.Body[n - 1 + k][k + 1] -= 0.5 * b2 * SinTetta[k + 1] / Math.Sqrt(1 - Math.Pow(_sArray[k + 1], 2));
                    _matrixA.Body[n - 1 + k][n + k] -= 0.5 * b2 * CosTetta[k] / Math.Sqrt(1 - Math.Pow(_sArray[k], 2));
                    _matrixA.Body[n - 1 + k][n + k + 1] -= 0.5 * b2 * CosTetta[k + 1] / Math.Sqrt(1 - Math.Pow(_sArray[k + 1], 2));
                }
                //Direct filling 1,2n
                for (int l = 0; l < n; l++)
                {
                    //Last Rows
                    _matrixA.Body[2 * n - 2][l] += CosTetta[l] / n;
                    _matrixA.Body[2 * n - 2][l + n] += -SinTetta[l] / n;
                    _matrixA.Body[2 * n - 1][l] += SinTetta[l] / n;
                    _matrixA.Body[2 * n - 1][l + n] += CosTetta[l] / n;


                    for (int m = 0; m < n - 1; m++)
                    {
                        //(1, n-1)
                        _matrixA.Body[m][l] += SinTetta[l] / (n * (_sArray[l] - _tArray[m])) 
                            + K11(_sArray[l], _tArray[m]) * CosTetta[l] / n
                            + K12(_sArray[l], _tArray[m]) * SinTetta[l] / n;
                        _matrixA.Body[m][l + n] += CosTetta[l] / (n * (_sArray[l] - _tArray[m]))
                            - K11(_sArray[l], _tArray[m]) * SinTetta[l] / n
                            + K12(_sArray[l], _tArray[m]) * CosTetta[l] / n;
                        //(n, 2n-2)
                        _matrixA.Body[m + n - 1][l] += CosTetta[l] / (n * (_sArray[l] - _tArray[m]))
                            + K21(_sArray[l], _tArray[m]) * CosTetta[l] / n
                            + K22(_sArray[l], _tArray[m]) * SinTetta[l] / n;
                        _matrixA.Body[m + n - 1][l + n] += -SinTetta[l] / (n * (_sArray[l] - _tArray[m]))
                            - K21(_sArray[l], _tArray[m]) * SinTetta[l] / n
                            + K22(_sArray[l], _tArray[m]) * CosTetta[l] / n;
                    }
                }
                #endregion
                //printer(p11 + Environment.NewLine);
                p11array[i] = p11;
                detarray[i] = _matrixA.Determinant;
                //printer(detarray[i] + Environment.NewLine);

                // when det < 0 stop
                if (detarray[i] < 0.0f)
                {
                    if (div > 0.00001)
                    {
                        dec = p11 - div;
                        div *= 0.1;
                        Calculate(n, h, a, nu1, nu2, mu1, mu2);
                        break;
                    }
                    //printer("*---------------------------*" + Environment.NewLine);
                    printer("h = " + h + Environment.NewLine);//p11 determ
                    printer($"Critical P {p11array[i - 1]}  {Environment.NewLine}");/*{detarray[i - 1]}*/
                    //printer($"{p11array[i]} {detarray[i]}{Environment.NewLine}");
                    //for (int last = i; last < 101; last++)
                       // backgroundWorker.ReportProgress(last);
                    //printer("*---------------------------*" + Environment.NewLine);
                    div = 0.01; dec = 0;
                    Printrepres(_matrixA.AbstractResolutions);
                    double[] answ = new double[23];
                    double[] tmp = _matrixA.AbstractResolutions;
                    if (mu1 < mu2)
                    {
                        for (int v = 0; v < 23; v++)
                        {
                            answ[v] = tmp[52 + v];
                        }
                    }
                    else
                    {
                        for (int v = 0; v < 23; v++)
                        {
                            answ[v] = tmp[97 - v];
                        }
                    }
                    for (int v = 0; v < 23; v++)
                    {
                        globalSolutione[23 + v] = answ[22 - v];
                        globalSolutione[v] = answ[v];
                    }
                    Printrepres(globalSolutione);
                    return globalSolutione;
                }
                backgroundWorker.ReportProgress((int)(i / dec / 2));
            }
            return globalSolutione;
        }

        private void Printrepres(double[] matrixAAbstractResolutions)
        {
            printer($"answrs = {Environment.NewLine}");
            
            for (int i = 0; i < 23; i++)
                printer($"{matrixAAbstractResolutions[i]},{Environment.NewLine}");
        }
    }
}
