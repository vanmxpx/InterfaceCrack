using System;

namespace InterfaceCrack
{
    public class Matrix
    {
        public int Size { get; }

        public double[][] Body { get; }

        public double[] Answer { get; }

        private bool issue = false;
        public double[] AbstractResolutions
        {
            get
            {
                if (issue) return Answer;
                 Answer[Size - 1] = 1;
                
                return ResolveLinearEquation();
            }
        }

        public double Determinant => CalculateDeterminant();

        public Matrix(int rows)
        {
            Size = rows;
            Answer = new double[rows];
            Body = new double[rows][];
            Null();
        }

        private double CalculateDeterminant()
        {
            // input checks
            double[] tmp = new double[Size];

            // pivoting
            for (int col = 0; col + 1 < Size; col++)
                if (Body[col][col] == 0)
                    // check for zero coefficients
                {
                    // find non-zero coefficient
                    int swapRow = col + 1;
                    for (; swapRow < Size; swapRow++) if (Body[swapRow][col] != 0) break;

                    if (Body[swapRow][col] != 0) //if found a non-zero coefficient
                    {
                        // swap it with the above
                        tmp = Body[swapRow];
                        Body[swapRow] = Body[col];
                        Body[col] = tmp;
                    }
                    else return 0; //else the matrix has no unique solution
                }
            // elimination
            double det = 1;
            for (int sourceRow = 0; sourceRow + 1 < Size; sourceRow++)
            {
                for (int destRow = sourceRow + 1; destRow < Size; destRow++)
                {
                    double df = Body[sourceRow][sourceRow];
                    double sf = Body[destRow][sourceRow];
                    for (int i = 0; i < Size; i++)
                        Body[destRow][i] -= Body[sourceRow][i] * sf / df;
                }
                det *= Body[sourceRow][sourceRow];
            }
            return det * Body[Size - 1][Size - 1];
        }

        private double[] ResolveLinearEquation()
        {
            issue = true;
            // back-insertion
            for (int row = Size - 1; row >= 0; row--)
            {
                double f = Body[row][row];
                if (f == 0) return new double[0];

                for (int i = 0; i < Size; i++)
                    Body[row][i] /= f;
                Answer[row] /= f;
                for (int destRow = 0; destRow < row; destRow++)
                {
                    Answer[destRow] -= Body[destRow][row] * Answer[row];
                    Body[destRow][row] = 0;
                }
            }

            return Answer;
        }

        public void Print(Action<string> printer, bool withSolution)
        {
            for (int i = 0; i < Size; i++)
            {
                for (int k = 0; k < Size; k++)
                {
                    printer(Body[i][k] + " ");
                }
                if (withSolution)
                    printer(Answer[i].ToString());
                printer(Environment.NewLine);
            }
            printer(Environment.NewLine);
        }

        public void Null()
        {
                issue = false;
            for (int i = 0; i < Size; i++)
                Body[i] = new double[Size];
        }
    }
}
