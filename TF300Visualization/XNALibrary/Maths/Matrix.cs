using System;
using System.Collections.Generic;
using System.Text;
using System.Collections;

namespace XNAHelper.Maths
{
    public class Matrix
    {

        /// <summary>
        /// Contains the rows of the matrix as elements, which
        /// are ArrayLists as well.
        /// </summary>
        private ArrayList Values;

        /// <summary>
        /// Number of rows of the matrix.
        /// </summary>
        public int RowCount
        {
            get { return rowCount; }
        }

        /// <summary>
        /// Number of columns of the matrix.
        /// </summary>
        public int ColumnCount
        {
            get { return columnCount; }
        }

        /// <summary>
        /// Number of rows of the matrix.
        /// </summary>
        private int rowCount;

        /// <summary>
        /// Number of columns of the matrix.
        /// </summary>
        private int columnCount;
        
        #region Constructors

        /// <summary>
        /// Inits empty matrix 
        /// </summary>
        public Matrix()
        {
            Values = new ArrayList();
            rowCount = 0;
            columnCount = 0;
        }

        /// <summary>
        /// Creates m by n matrix filled with zeros; same as Zeros(m, n).
        /// </summary>
        /// <param name="m">Number of rows</param>
        /// <param name="n">Number of columns</param>
        public Matrix(int m, int n)
        {
            rowCount = m;
            columnCount = n;

            Values = new ArrayList(m);

            for (int i = 0; i < m; i++)
            {
                Values.Add(new ArrayList(n));

                for (int j = 0; j < n; j++)
                {
                    ((ArrayList)Values[i]).Add(Complex.Zero);
                }
            }
        }

        /// <summary>
        /// Inits square matrix
        /// </summary>
        /// <param name="n"></param>
        public Matrix(int n)
        {
            rowCount = n;
            columnCount = n;

            Values = new ArrayList(n);

            for (int i = 0; i < n; i++)
            {
                Values.Add(new ArrayList(n));

                for (int j = 0; j < n; j++)
                {
                    ((ArrayList)Values[i]).Add(Complex.Zero);
                }
            }
        }

        /// <summary>
        /// Creates one by one matrix containing x
        /// </summary>
        /// <param name="x"></param>
        public Matrix(Complex x)
        {
            rowCount = 1;
            columnCount = 1;

            Values = new ArrayList(1);

            Values.Add(new ArrayList(1));

            ((ArrayList)Values[0]).Add(x);
        }

        /// <summary>
        /// Creates matrix from 2-d Complex array.
        /// </summary>
        /// <param name="values"></param>
        public Matrix(Complex[,] values)
        {
            if (values == null)
            {
                Values = new ArrayList();
                columnCount = 0;
                rowCount = 0;
            }

            rowCount = (int)values.GetLongLength(0);
            columnCount = (int)values.GetLongLength(1);

            Values = new ArrayList(rowCount);

            for (int i = 0; i < rowCount; i++)
            {
                Values.Add(new ArrayList(columnCount));

                for (int j = 0; j < columnCount; j++)
                {
                    ((ArrayList)Values[i]).Add(values[i, j]);
                }
            }
        }

        /// <summary>
        /// Creates column vector from Complex array.
        /// </summary>
        /// <param name="values"></param>
        public Matrix(Complex[] values)
        {
            if (values == null)
            {
                Values = new ArrayList();
                columnCount = 0;
                rowCount = 0;
            }

            rowCount = values.Length;
            columnCount = 1;

            Values = new ArrayList(rowCount);

            for (int i = 0; i < rowCount; i++)
            {
                Values.Add(new ArrayList(1));

                ((ArrayList)Values[i]).Add(values[i]);
            }
        }

        /// <summary>
        /// Creates one by one matrix containing x
        /// </summary>
        /// <param name="x"></param>
        public Matrix(double x)
        {
            rowCount = 1;
            columnCount = 1;

            Values = new ArrayList(1);

            Values.Add(new ArrayList(1));

            ((ArrayList)Values[0]).Add(new Complex(x));
        }

        /// <summary>
        /// Creates matrix from 2-d double array.
        /// </summary>
        /// <param name="values"></param>
        public Matrix(double[,] values)
        {
            if (values == null)
            {
                Values = new ArrayList();
                columnCount = 0;
                rowCount = 0;
            }

            rowCount = (int)values.GetLongLength(0);
            columnCount = (int)values.GetLongLength(1);

            Values = new ArrayList(rowCount);

            for (int i = 0; i < rowCount; i++)
            {
                Values.Add(new ArrayList(columnCount));

                for (int j = 0; j < columnCount; j++)
                {
                    ((ArrayList)Values[i]).Add(new Complex(values[i, j]));
                }
            }
        }

        /// <summary>
        /// Creates column vector from double array.
        /// </summary>
        /// <param name="values"></param>
        public Matrix(double[] values)
        {
            if (values == null)
            {
                Values = new ArrayList();
                columnCount = 0;
                rowCount = 0;
            }

            rowCount = values.Length;
            columnCount = 1;

            Values = new ArrayList(rowCount);

            for (int i = 0; i < rowCount; i++)
            {
                Values.Add(new ArrayList(1));

                ((ArrayList)Values[i]).Add(new Complex(values[i]));
            }
        }

        /// <summary>
        /// Creates real matrix from string, e.g. "1,0;0,1" gives the 2 by 2 identity matrix.
        /// Not fast, but easy to use, if matrices are to be entered by hand or read from text files.
        /// </summary>
        /// <param name="matrix">Matrix coded as string. Lines are separated by a semicolon, column elements by a comma.</param>
        public Matrix(string matrix_string)
        {
            // remove spaces
            matrix_string = matrix_string.Replace(" ", "");

            // split string into rows, use ';' as separator
            string[] rows = matrix_string.Split(new char[] { ';' });

            // init Values, RowCount, ColumnCount
            rowCount = rows.Length;
            Values = new ArrayList(rowCount);
            columnCount = 0;

            for (int i = 0; i < rowCount; i++)
                Values.Add(new ArrayList());

            string[] curcol;

            //try
            {
                for (int i = 1; i <= rowCount; i++)
                {
                    curcol = rows[i - 1].Split(new char[] { ',' });

                    for (int j = 1; j <= curcol.Length; j++)
                    {
                        this[i, j] = new Complex(Convert.ToDouble(curcol[j - 1]));
                    }
                }
            }

        }

        #endregion

        #region Static func

        /// <summary>
        /// Retrieves the j-th canoncical basis vector of the IR^n.
        /// </summary>
        /// <param name="n">Dimension of the basis.</param>
        /// <param name="j">Index of canonical basis vector to be retrieved.</param>
        /// <returns></returns>
        public static Matrix E(int n, int j)
        {
            Matrix e = Zeros(n, 1);
            e[j] = Complex.One;

            return e;
        }

        /// <summary>
        /// Returns 1 if i = j, and 0 else.
        /// </summary>
        /// <param name="i"></param>
        /// <param name="j"></param>
        /// <returns></returns>
        public static Complex KroneckerDelta(int i, int j)
        {
            return new Complex(Math.Min(Math.Abs(i - j), 1));
        }

        /// <summary>
        /// Creates m by n chessboard matrix with interchang韓g ones and zeros.
        /// 
        /// </summary>
        /// <param name="m">Number of rows.</param>
        /// <param name="n">Number of columns.</param>
        /// <param name="even">Indicates, if matrix entry (1,1) equals zero.</param>
        /// <returns></returns>
        public static Matrix ChessboardMatrix(int m, int n, bool even)
        {
            Matrix M = new Matrix(m, n);
            
            if (even)
                for (int i = 1; i <= m; i++)
                    for (int j = 1; j <= n; j++)
                        M[i, j] = KroneckerDelta((i + j) % 2, 0);
            else
                for (int i = 1; i <= m; i++)
                    for (int j = 1; j <= n; j++)
                        M[i, j] = KroneckerDelta((i + j) % 2, 1);

            return M;
        }

        /// <summary>
        /// Creates m by n chessboard matrix with interchang韓g ones and zeros.
        /// 
        /// </summary>        
        /// <param name="n">Number of columns.</param>
        /// <param name="even">Indicates, if matrix entry (1,1) equals zero.</param>
        /// <returns></returns>
        public static Matrix ChessboardMatrix(int n, bool even)
        {
            Matrix M = new Matrix(n);

            if (even)
                for (int i = 1; i <= n; i++)
                    for (int j = 1; j <= n; j++)
                        M[i, j] = KroneckerDelta((i + j) % 2, 0);
            else
                for (int i = 1; i <= n; i++)
                    for (int j = 1; j <= n; j++)
                        M[i, j] = KroneckerDelta((i + j) % 2, 1);

            return M;
        }

        /// <summary>
        /// Creates m by n matrix filled with zeros.
        /// </summary>
        /// <param name="m">Number of rows.</param>
        /// <param name="n">Number of columns.</param>
        /// <returns>m by n matrix filled with zeros.</returns>
        public static Matrix Zeros(int m, int n)
        {
            return new Matrix(m, n);
        }

        /// <summary>
        /// Creates n by n matrix filled with zeros.
        /// </summary>       
        /// <param name="n">Number of rows and columns, resp.</param>
        /// <returns>n by n matrix filled with zeros.</returns>
        public static Matrix Zeros(int n)
        {
            return new Matrix(n);
        }

        /// <summary>
        /// Creates m by n matrix filled with ones.
        /// </summary>
        /// <param name="m">Number of rows.</param>
        /// <param name="n">Number of columns.</param>
        /// <returns>m by n matrix filled with ones.</returns>        
        public static Matrix Ones(int m, int n)
        {
            Matrix M = new Matrix(m, n);

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    ((ArrayList)M.Values[i])[j] = Complex.One;
                }
            }

            return M;
        }

        /// <summary>
        /// Creates n by n matrix filled with ones.
        /// </summary>        
        /// <param name="n">Number of columns.</param>
        /// <returns>n by n matrix filled with ones.</returns>        
        public static Matrix Ones(int n)
        {
            Matrix M = new Matrix(n);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    ((ArrayList)M.Values[i])[j] = Complex.One;
                }
            }

            return M;
        }

        /// <summary>
        /// Creates n by n identity matrix.
        /// </summary>
        /// <param name="n">Number of rows and columns respectively.</param>
        /// <returns>n by n identity matrix.</returns>
        public static Matrix Identity(int n)
        {
            return Diag(Ones(n, 1));
        }

        /// <summary>
        /// Creates teh n by n identity matrix.
        /// </summary>
        /// <param name="n">Number of rows and columns, resp.</param>
        /// <returns></returns>
        public static Matrix Eye(int n)
        {
            return Identity(n);
        }

        /// <summary>
        /// Vertically concats matrices A and B, which do not have to be of the same height. 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns>Matrix [A|B]</returns>
        public static Matrix VerticalConcat(Matrix A, Matrix B)
        {

            Matrix C = A.Column(1);
            
            for (int j = 2; j <= A.ColumnCount; j++)
            {
                C.InsertColumn(A.Column(j), j);
            }


            for (int j = 1; j <= B.ColumnCount; j++)
            {
                C.InsertColumn(B.Column(j), C.ColumnCount + 1);
            }

            return C;
        }

        public static Matrix VerticalConcat(Matrix[] A)
        {
            if (A == null)
                throw new ArgumentNullException();
            else if (A.Length == 1)
                return A[0];
            else
            {
                Matrix C = VerticalConcat(A[0], A[1]);

                for (int i = 2; i < A.Length; i++)
                {
                    C = VerticalConcat(C, A[i]);
                }

                return C;
            }
        }

        public static Matrix HorizontalConcat(Matrix A, Matrix B)
        {
            Matrix C = A.Row(1);

          
            for (int i = 2; i <= A.RowCount; i++)
            {
                C.InsertRow(A.Row(i), i);
            }


            for (int i = 1; i <= B.RowCount; i++)
            {
                C.InsertRow(B.Row(i), C.RowCount + 1);
            }

            return C;
        }

        public static Matrix HorizontalConcat(Matrix[] A)
        {
            if (A == null)
                throw new ArgumentNullException();
            else if (A.Length == 1)
                return A[0];
            else
            {
                Matrix C = HorizontalConcat(A[0], A[1]);

                for (int i = 2; i < A.Length; i++)
                {
                    C = HorizontalConcat(C, A[i]);
                }

                return C;
            }
        }


        /// <summary>
        /// Generates diagonal matrix 计算对角矩阵 参见http://zh.wikipedia.org/zh-cn/%E5%B0%8D%E8%A7%92%E7%9F%A9%E9%99%A3  
        /// </summary>
        /// <param name="diag_vector">column vector containing the diag elements</param>
        /// <returns></returns>
        public static Matrix Diag(Matrix diag_vector)
        {
            int dim = diag_vector.VectorLength();
            
            if(dim == 0)
                throw new ArgumentException("diag_vector must be 1xN or Nx1");

            Matrix M = new Matrix(dim, dim);

            for (int i = 1; i <= dim; i++)
            {
                M[i, i] = diag_vector[i];
            }

            return M;

        }

        /// <summary>
        /// Generates diagonal matrix
        /// </summary>
        /// <param name="diag_vector">column vector containing the diag elements</param>
        /// <returns></returns>
        public static Matrix Diag(Matrix diag_vector, int offset)
        {
            int dim = diag_vector.VectorLength();
            

            if(dim == 0)
                throw new ArgumentException("diag_vector must be 1xN or Nx1.");

            //if (Math.Abs(offset) >= dim)
            //    throw new ArgumentException("Absolute value of offset must be less than length of diag_vector.");

            Matrix M = new Matrix(dim + Math.Abs(offset), dim + Math.Abs(offset));
            dim = M.RowCount;

            if (offset >= 0)
            {
                for (int i = 1; i <= dim - offset; i++)
                {
                    M[i, i + offset] = diag_vector[i];
                }
            }
            else
            {
                for (int i = 1; i <= dim + offset; i++)
                {
                    M[i - offset, i] = diag_vector[i];
                }
            }

            return M;

        }

        /// <summary>
        /// Generates tri-diagonal square matrix with constant values on main
        /// and secondary diagonals.
        /// </summary>
        /// <param name="l">Value of lower secondary diagonal.</param>
        /// <param name="d">Value of main diagonal.</param>
        /// <param name="u">Value of upper secondary diagonal.</param>
        /// <param name="n">Dimension of the output matrix.</param>
        /// <returns>nxn tri-diagonal matrix.</returns>
        public static Matrix TriDiag(Complex l, Complex d, Complex u, int n)
        {
            if (n <= 1)
                throw new ArgumentException("Matrix dimension must greater than one.");

            return Diag(l * Ones(n - 1, 1), -1) + Diag(d * Ones(n, 1)) + Diag(u * Ones(n - 1, 1), 1);
        }

        /// <summary>
        /// Generates tri-diagonal square matrix with overloaded vectors
        /// as main and secondary diagonals. The dimension of the output
        /// matrix is determined by the length of d.
        /// </summary>
        /// <param name="l">Lower secondary diagonal vector.</param>
        /// <param name="d">Main diagonal vector.</param>
        /// <param name="u">Upper secondary diagonal vector.</param>
        /// <returns></returns>
        public static Matrix TriDiag(Matrix l, Matrix d, Matrix u)
        {
            int sizeL = l.VectorLength();
            int sizeD = d.VectorLength();
            int sizeU = u.VectorLength();

            if (sizeL * sizeD * sizeU == 0)
                throw new ArgumentException("At least one of the paramter matrices is not a vector.");

            if (sizeL != sizeU)
                throw new ArgumentException("Lower and upper secondary diagonal must have the same length.");

            if (sizeL + 1 != sizeD)
                throw new ArgumentException("Main diagonal must have exactly one element more than the secondary diagonals.");

            return Diag(l, -1) + Diag(d) + Diag(u, 1);
        }



        /// <summary>
        /// Implements the dot product of two vectors.
        /// </summary>
        /// <param name="v">Row or column vector.</param>
        /// <param name="w">Row or column vector.</param>
        /// <returns>Dot product.</returns>
        public static Complex Dot(Matrix v, Matrix w)
        {
            int m = v.VectorLength();
            int n = w.VectorLength();

            if (m == 0 || n == 0)
                throw new ArgumentException("Arguments need to be vectors.");
            else if (m != n)
                throw new ArgumentException("Vectors must be of the same length.");

            Complex buf = Complex.Zero;

            for (int i = 1; i <= m; i++)
            {
                buf += v[i] * w[i];
            }

            return buf;
        }

        /// <summary>
        /// Calcs the n-th Fibonacci-number in O(n)
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static Complex Fib(int n)
        {
            Matrix M = Ones(2, 2);
            M[2, 2] = Complex.Zero;

            return (M ^ (n - 1))[1, 1];
        }

        /// <summary>
        /// Creates n by n matrix filled with random values in [0,1];
        /// all entries on the main diagonal are zero.
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static Matrix RandomGraph(int n)
        {
            Matrix buf = Random(n, n);

            buf -= Diag(buf.DiagVector());

            return buf;
        }

        /// <summary>
        /// Creates n by n matrix filled with random values in [0,1];
        /// all entries on the main diagonal are zero.
        /// A specified random percentage of edges has weight positive infinity.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="p">Defines probability for an edge being less than +infty. Should be in [0,1],
        /// p = 1 gives complete directed graph; p = 0 gives no edges.</param>
        /// <returns></returns>
        public static Matrix RandomGraph(int n, double p)
        {
            Matrix buf = new Matrix(n);

            Random r = new Random();

            for (int i = 1; i <= n; i++)
                for (int j = 1; j <= n; j++)
                    if (i == j) buf[i, j] = Complex.Zero;
                    else if (r.NextDouble() < p) buf[i, j] = new Complex(r.NextDouble());
                    else buf[i, j] = new Complex(double.PositiveInfinity);

            return buf;
        }

        /// <summary>
        /// Creates m by n matrix filled with random values in [0,1].
        /// </summary>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static Matrix Random(int m, int n)
        {
            Matrix M = new Matrix(m, n);
            Random r = new Random();

            for (int i = 1; i <= m; i++)
            {
                for (int j = 1; j <= n; j++)
                {
                    M[i, j] = new Complex(r.NextDouble());
                }
            }

            return M;
        }

        /// <summary>
        /// Creates n by n matrix filled with random values in [0,1].
        /// </summary>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static Matrix Random(int n)
        {
            Matrix M = new Matrix(n);
            Random r = new Random();

            for (int i = 1; i <= n; i++)
            {
                for (int j = 1; j <= n; j++)
                {
                    M[i, j] = new Complex(r.NextDouble());
                }
            }

            return M;
        }

        /// <summary>
        /// Creates n by n matrix filled with random values in {lo,...,hi-1}.
        /// </summary>
        ///<param name="lo">Inclusive lower bound.</param>
        /// <param name="hi">Exclusive upper bound</param>
        /// <param name="n">Number of rows and columns each.</param>
        /// <returns></returns>
        public static Matrix Random(int n, int lo, int hi)
        {
            Matrix M = new Matrix(n);
            Random r = new Random();

            for (int i = 1; i <= n; i++)
            {
                for (int j = 1; j <= n; j++)
                {
                    M[i, j] = new Complex((double)r.Next(lo, hi));
                }
            }

            return M;
        }

        /// <summary>
        /// Creates m by n random zero one matrix with probability p for a one.
        /// </summary>
        /// <param name="m">Number of rows.</param>
        /// <param name="n">Number of columns.</param>
        /// <param name="p">Probability fro an entry to be one, expecting a value in [0,1].</param>
        /// <returns></returns>
        public static Matrix RandomZeroOne(int m, int n, double p)
        {
            Matrix M = new Matrix(m, n);
            Random r = new Random();

            for (int i = 1; i <= m; i++)
                for (int j = 1; j <= n; j++)
                    if (r.NextDouble() <= p) M[i, j] = Complex.One;

            return M;
        }

        /// <summary>
        /// Creates n by n random zero one matrix with probability p for a one.
        /// </summary>        
        /// <param name="n">Number of rows and columns, resp.</param>
        /// <param name="p">Probability fro an entry to be one, expecting a value in [0,1].</param>
        /// <returns></returns>
        public static Matrix RandomZeroOne(int n, double p)
        {
            Matrix M = new Matrix(n, n);
            Random r = new Random();

            for (int i = 1; i <= n; i++)
                for (int j = 1; j <= n; j++)
                    if (r.NextDouble() <= p) M[i, j] = Complex.One;

            return M;
        }

        /// <summary>
        /// Creates m by n matrix filled with random values in {lo,...,hi-1}.
        /// </summary>
        ///<param name="lo">Inclusive lower bound.</param>
        /// <param name="hi">Exclusive upper bound</param>
        /// <param name="m">Number of rows.</param>
        /// <param name="n">Number of columns.</param>
        /// <returns></returns>
        public static Matrix Random(int m, int n, int lo, int hi)
        {
            Matrix M = new Matrix(m, n);
            Random r = new Random();

            for (int i = 1; i <= m; i++)
            {
                for (int j = 1; j <= n; j++)
                {
                    M[i, j] = new Complex((double)r.Next(lo, hi));
                }
            }

            return M;
        }

        public static Matrix Vandermonde(Complex[] x)
        {
            if (x == null || x.Length < 1)
                throw new ArgumentNullException();

            int n = x.Length - 1;

            Matrix V = new Matrix(n + 1);

            for (int i = 0; i <= n; i++)
                for (int p = 0; p <= n; p++)
                    V[i + 1, p + 1] = Complex.Pow(x[i], p);

            return V;
        }

        /// <summary>
        /// Computes all shortest distance between any vertices in a given graph.
        /// </summary>
        /// <param name="adjacence_matrix">Square adjacence matrix. The main diagonal
        /// is expected to consist of zeros, any non-existing edges should be marked
        /// positive infinity.</param>
        /// <returns>Two matrices D and P, where D[u,v] holds the distance of the shortest
        /// path between u and v, and P[u,v] holds the shortcut vertex on the way from
        /// u to v.</returns>
        public static Matrix[] Floyd(Matrix adjacence_matrix)
        {
            if (!adjacence_matrix.IsSquare())
                throw new ArgumentException("Expected square matrix.");
            else if (!adjacence_matrix.IsReal())
                throw new ArgumentException("Adjacence matrices are expected to be real.");

            int n = adjacence_matrix.RowCount;

            Matrix D = adjacence_matrix.Clone(); // distance matrix
            Matrix P = new Matrix(n);

            double buf;

            for (int k = 1; k <= n; k++)
                for (int i = 1; i <= n; i++)
                    for (int j = 1; j <= n; j++)
                    {
                        buf = D[i, k].Re + D[k, j].Re;
                        if (buf < D[i, j].Re)
                        {
                            D[i, j].Re = buf;
                            P[i, j].Re = k;
                        }
                    }

            return new Matrix[] { D, P };
        }

        /// <summary>
        /// Returns the shortest path between two given vertices i and j as
        /// int array.
        /// </summary>
        /// <param name="P">Path matrix as returned from Floyd().</param>
        /// <param name="i">One-based index of start vertex.</param>
        /// <param name="j">One-based index of end vertex.</param>
        /// <returns></returns>
        public static ArrayList FloydPath(Matrix P, int i, int j)
        {
            if (!P.IsSquare())
                throw new ArgumentException("Path matrix must be square.");
            else if (!P.IsReal())
                throw new ArgumentException("Adjacence matrices are expected to be real.");

            ArrayList path = new ArrayList();
            path.Add(i);

            //int borderliner = 0;
            //int n = P.Size()[0] + 1; // shortest path cannot have more than n vertices! 

            while (P[i, j] != 0)
            {
                i = Convert.ToInt32(P[i, j]);
                path.Add(i);

                //borderliner++;

                //if (borderliner == n)
                //    throw new FormatException("P was not a Floyd path matrix.");
            }

            path.Add(j);

            return path;
        }

        /// <summary>
        /// Performs depth-first search for a graph given by its adjacence matrix.
        /// </summary>
        /// <param name="adjacence_matrix">A[i,j] = 0 or +infty, if there is no edge from i to j; any non-zero value otherwise.</param>
        /// <param name="root">The vertex to begin the search.</param>
        /// <returns>Adjacence matrix of the computed spanning tree.</returns>
        public static Matrix DFS(Matrix adjacence_matrix, int root)
        {
            if (!adjacence_matrix.IsSquare())
                throw new ArgumentException("Adjacence matrices are expected to be square.");
            else if (!adjacence_matrix.IsReal())
                throw new ArgumentException("Adjacence matrices are expected to be real.");


            int n = adjacence_matrix.RowCount;


            if (root < 1 || root > n)
                throw new ArgumentException("Root must be a vertex of the graph, e.i. in {1, ..., n}.");

            Matrix spanTree = new Matrix(n);

            bool[] marked = new bool[n + 1];

            Stack todo = new Stack();
            todo.Push(root);
            marked[root] = true;

            // adajacence lists for each vertex
            ArrayList[] A = new ArrayList[n + 1];

            for (int i = 1; i <= n; i++)
            {
                A[i] = new ArrayList();

                for (int j = 1; j <= n; j++)
                    if (adjacence_matrix[i, j].Re != 0 && adjacence_matrix[i, j].Im != double.PositiveInfinity)
                        A[i].Add(j);
            }

            int v, w;

            while (todo.Count > 0)
            {
                v = (int)todo.Peek();

                if (A[v].Count > 0)
                {
                    w = (int)A[v][0];

                    if (!marked[w])
                    {
                        marked[w] = true; // mark w
                        spanTree[v, w].Re = 1; // mark vw
                        todo.Push(w); // one more to search
                    }

                    A[v].RemoveAt(0);
                }
                else
                    todo.Pop();
            }

            return spanTree;
        }

        /// <summary>
        /// Performs broad-first search for a graph given by its adjacence matrix.
        /// </summary>
        /// <param name="adjacence_matrix">A[i,j] = 0 or +infty, if there is no edge from i to j; any non-zero value otherwise.</param>
        /// <param name="root">The vertex to begin the search.</param>
        /// <returns>Adjacence matrix of the computed spanning tree.</returns>
        public static Matrix BFS(Matrix adjacence_matrix, int root)
        {
            if (!adjacence_matrix.IsSquare())
                throw new ArgumentException("Adjacence matrices are expected to be square.");
            else if (!adjacence_matrix.IsReal())
                throw new ArgumentException("Adjacence matrices are expected to be real.");

            int n = adjacence_matrix.RowCount;


            if (root < 1 || root > n)
                throw new ArgumentException("Root must be a vertex of the graph, e.i. in {1, ..., n}.");

            Matrix spanTree = new Matrix(n);

            bool[] marked = new bool[n + 1];

            Queue todo = new Queue();
            todo.Enqueue(root);
            marked[root] = true;

            // adajacence lists for each vertex
            ArrayList[] A = new ArrayList[n + 1];

            for (int i = 1; i <= n; i++)
            {
                A[i] = new ArrayList();

                for (int j = 1; j <= n; j++)
                    if (adjacence_matrix[i, j].Re != 0 && adjacence_matrix[i, j].Re != double.PositiveInfinity)
                        A[i].Add(j);
            }

            int v, w;

            while (todo.Count > 0)
            {
                v = (int)todo.Peek();

                if (A[v].Count > 0)
                {
                    w = (int)A[v][0];

                    if (!marked[w])
                    {
                        marked[w] = true; // mark w
                        spanTree[v, w].Re = 1; // mark vw
                        todo.Enqueue(w); // one more to search
                    }

                    A[v].RemoveAt(0);
                }
                else
                    todo.Dequeue();
            }

            return spanTree;
        }

        /// <summary>
        /// Creates a random matrix filled with zeros and ones.
        /// </summary>
        /// <param name="m">Number of rows.</param>
        /// <param name="n">Number of columns.</param>
        /// <param name="p">Probability of each entry being 1.</param>
        /// <returns></returns>
        public static Matrix ZeroOneRandom(int m, int n, double p)
        {
            Random r = new Random();

            Matrix buf = Zeros(m, n);

            for (int i = 1; i <= m; i++)
            {
                for (int j = 1; j <= n; j++)
                {
                    if (r.NextDouble() <= p)
                        buf[i, j] = Complex.One;
                }
            }

            return buf;
        }

        /// <summary>
        /// Creates a random matrix filled with zeros and ones.
        /// </summary>        
        /// <param name="n">Number of rows and columns.</param>
        /// <param name="p">Probability of each entry being 1.</param>
        /// <returns></returns>
        public static Matrix ZeroOneRandom(int n, double p)
        {
            Random r = new Random();

            Matrix buf = Zeros(n);

            for (int i = 1; i <= n; i++)
            {
                for (int j = 1; j <= n; j++)
                {
                    if (r.NextDouble() <= p)
                        buf[i, j] = Complex.One;
                }
            }

            return buf;
        }

        /// <summary>
        /// Computes the Householder vector.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private static Matrix[] HouseholderVector(Matrix x)
        {
            //throw new NotImplementedException("Supposingly buggy!");

            //if (!x.IsReal())
            //    throw new ArgumentException("Cannot compute housholder vector of non-real vector.");

            int n = x.VectorLength();

            if (n == 0)
                throw new InvalidOperationException("Expected vector as argument.");

            Matrix y = x / x.Norm();
            Matrix buf = y.Extract(2, n, 1, 1);
            Complex s = Dot(buf, buf);

            Matrix v = Zeros(n, 1);
            v[1] = Complex.One;

            v.Insert(2, 1, buf);

            double beta = 0;

            if (s != 0)
            {
                Complex mu = Complex.Sqrt(y[1] * y[1] + s);
                if (y[1].Re <= 0)
                    v[1] = y[1] - mu;
                else
                    v[1] = -s / (y[1] + mu);

                beta = 2 * v[1].Re * v[1].Re / (s.Re + v[1].Re * v[1].Re);
                v = v / v[1];
            }

            return new Matrix[] { v, new Matrix(beta) };
            
        }

        /// <summary>
        /// Constructs block matrix [A, B; C, D].
        /// </summary>
        /// <param name="A">Upper left sub matrix.</param>
        /// <param name="B">Upper right sub matrix.</param>
        /// <param name="C">Lower left sub matrix.</param>
        /// <param name="D">Lower right sub matrix.</param>
        /// <returns></returns>
        public static Matrix BlockMatrix(Matrix A, Matrix B, Matrix C, Matrix D)
        {
           
            if (A.RowCount != B.RowCount || C.RowCount != D.RowCount
                || A.ColumnCount != C.ColumnCount || B.ColumnCount != D.ColumnCount)
                throw new ArgumentException("Matrix dimensions must agree.");

            Matrix R = new Matrix(A.RowCount + C.RowCount, A.ColumnCount + B.ColumnCount);

            for (int i = 1; i <= R.rowCount; i++)
                for (int j = 1; j <= R.columnCount; j++)
                    if (i <= A.RowCount)
                    {
                        if (j <= A.ColumnCount)
                            R[i, j] = A[i, j];
                        else
                            R[i, j] = B[i, j - A.ColumnCount];
                    }
                    else
                    {
                        if (j <= C.ColumnCount)
                            R[i, j] = C[i - A.RowCount, j];
                        else
                            R[i, j] = D[i - A.RowCount, j - C.ColumnCount];
                    }

            return R;
        }

        /// <summary>
        /// For this matrix A, this method solves Ax = b via LU factorization with
        /// column pivoting.
        /// </summary>
        /// <param name="b">Vector of appropriate length.</param>
        /// <remarks>Approximately n^3/3 + 2n^2 dot operations ~> O(n^3)</remarks>
        public static Matrix Solve(Matrix A, Matrix b)
        {
            Matrix A2 = A.Clone();
            Matrix b2 = b.Clone();


            if (!A2.IsSquare())
                throw new InvalidOperationException("Cannot uniquely solve non-square equation system.");

            int n = A2.RowCount;

            Matrix P = A2.LUSafe();

            // We know: PA = LU => [ Ax = b <=> P'LUx = b <=> L(Ux) = (Pb)] since P is orthogonal
            // set y := Ux, solve Ly = Pb by forward insertion
            // and Ux = y by backward insertion

            b2 = P * b2;
            
            // this solves Ly = Pb
            (A2.ExtractLowerTrapeze() - Diag(A2.DiagVector()) + Identity(n)).ForwardInsertion(b2);

            // this solves Ux = y
            (A2.ExtractUpperTrapeze()).BackwardInsertion(b2);

            return b2;

            
        }

        #endregion

        #region Dynamic funcs

        #region Matrix manipulations, extractions and decompositions

        /// <summary>
        /// Returns the matrix of the real parts of the entries of this matrix.
        /// </summary>
        /// <returns></returns>
        public Matrix Re()
        {
            Matrix M = new Matrix(rowCount, columnCount);

            for (int i = 1; i <= rowCount; i++)
                for (int j = 1; j <= columnCount; j++)
                    M[i, j] = new Complex(this[i, j].Re);

            return M;
        }

        /// <summary>
        /// Returns the matrix of the imaginary parts of the entries of this matrix.
        /// </summary>
        /// <returns></returns>
        public Matrix Im()
        {
            Matrix M = new Matrix(rowCount, columnCount);

            for (int i = 1; i <= rowCount; i++)
                for (int j = 1; j <= columnCount; j++)
                    M[i, j] = new Complex(this[i, j].Im);

            return M;
        }

        /// <summary>
        /// Performs Hessenberg-Householder reduction, where {H, Q}
        /// is returned, with H Hessenbergian, Q orthogonal and H = Q'AQ.
        /// </summary>
        /// <returns></returns>
        public Matrix[] HessenbergHouseholder()
        {
            //throw new NotImplementedException("Still buggy!");

            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot perform Hessenberg Householder decomposition of non-square matrix.");

            int n = rowCount;
            Matrix Q = Identity(n);
            Matrix H = this.Clone();
            Matrix I, N, R, P;
            Matrix[] vbeta = new Matrix[2];
            int m;

            // don't try to understand from the code alone.
            // this is pure magic to me - mathematics, reborn as code.
            for (int k = 1; k <= n - 2; k++)
            {
                vbeta = HouseholderVector(H.Extract(k + 1, n, k, k));
                I = Identity(k);
                N = Zeros(k, n - k);

                m = vbeta[0].VectorLength();
                R = Identity(m) - vbeta[1][1, 1] * vbeta[0] * vbeta[0].Transpose();

                H.Insert(k + 1, k, R * H.Extract(k + 1, n, k, n));
                H.Insert(1, k + 1, H.Extract(1, n, k + 1, n) * R);

                P = BlockMatrix(I, N, N.Transpose(), R);

                Q = Q * P;
            }

            return new Matrix[] { H, Q };

        }

        /// <summary>
        /// Extract sub matrix.
        /// </summary>
        /// <param name="i1">Start row.</param>
        /// <param name="i2">End row.</param>
        /// <param name="j1">Start column.</param>
        /// <param name="j2">End column.</param>
        /// <returns></returns>
        public Matrix Extract(int i1, int i2, int j1, int j2)
        {
            if (i2 < i1 || j2 < j1 || i1 <= 0 || j2 <= 0 || i2 > rowCount || j2 > columnCount)
                throw new ArgumentException("Index exceeds matrix dimension.");

            Matrix B = new Matrix(i2 - i1 + 1, j2 - j1 + 1);

            for (int i = i1; i <= i2; i++)
                for (int j = j1; j <= j2; j++)
                    B[i - i1 + 1, j - j1 + 1] = this[i, j];

            return B;
        }

        /// <summary>
        /// Extracts lower trapeze matrix of this matrix.
        /// </summary>
        /// <returns></returns>
        public Matrix ExtractLowerTrapeze()
        {
            Matrix buf = new Matrix(rowCount, columnCount);

            for (int i = 1; i <= rowCount; i++)
            {
                for (int j = 1; j <= i; j++)
                {
                    buf[i, j] = this[i, j];
                }
            }

            return buf;
        }

        /// <summary>
        /// Extracts upper trapeze matrix of this matrix.
        /// </summary>
        /// <returns></returns>
        public Matrix ExtractUpperTrapeze()
        {
            Matrix buf = new Matrix(rowCount, columnCount);

            for (int i = 1; i <= rowCount; i++)
            {
                for (int j = i; j <= columnCount; j++)
                {
                    buf[i, j] = this[i, j];
                }
            }

            return buf;
        }

        /// <summary>
        /// Splits matrix into its column vectors.
        /// </summary>
        /// <returns>Array of column vectors.</returns>
        public Matrix[] ColumnVectorize()
        {
            Matrix[] buf = new Matrix[columnCount];

            for (int j = 1; j <= buf.Length; j++)
            {
                buf[j] = this.Column(j);
            }

            return buf;
        }

        /// <summary>
        /// Splits matrix into its row vectors.
        /// </summary>
        /// <returns>Array of row vectors.</returns>
        public Matrix[] RowVectorize()
        {
            Matrix[] buf = new Matrix[rowCount];

            for (int i = 1; i <= buf.Length; i++)
            {
                buf[i] = this.Row(i);
            }

            return buf;
        }

        /// <summary>
        /// Flips matrix vertically.
        /// </summary>
        public void VerticalFlip()
        {
            Values.Reverse();
        }

        /// <summary>
        /// Flips matrix horizontally.
        /// </summary>
        public void HorizontalFlip()
        {
            for (int i = 0; i < rowCount; i++)
            {
                ((ArrayList)Values[i]).Reverse();
            }
        }

        /// <summary>
        /// Swaps columns at specified indices. The latter do not have to be ordered.
        /// When equal, nothing is done.
        /// </summary>
        /// <param name="j1">One-based index of first col.</param>
        /// <param name="j2">One-based index of second col.</param>       
        public void SwapColumns(int j1, int j2)
        {
            if (j1 <= 0 || j1 > columnCount || j2 <= 0 || j2 > columnCount)
                throw new ArgumentException("Indices must be positive and <= number of cols.");

            if (j1 == j2)
                return;

            // ArrayList indices are zero-based
            j1--;
            j2--;
            object buf;

            for (int i = 0; i < rowCount; i++)
            {
                buf = ((ArrayList)Values[i])[j1];
                ((ArrayList)Values[i])[j1] = ((ArrayList)Values[i])[j2];
                ((ArrayList)Values[i])[j2] = buf;
            }
        }

        /// <summary>
        /// Swaps rows at specified indices. The latter do not have to be ordered.
        /// When equal, nothing is done.
        /// </summary>
        /// <param name="i1">One-based index of first row.</param>
        /// <param name="i2">One-based index of second row.</param>        
        public void SwapRows(int i1, int i2)
        {
            if (i1 <= 0 || i1 > rowCount || i2 <= 0 || i2 > rowCount)
                throw new ArgumentException("Indices must be positive and <= number of rows.");

            if (i1 == i2)
                return;

            ArrayList buf = (ArrayList)Values[--i1];
            Values[i1] = Values[--i2];
            Values[i2] = buf;
        }

        /// <summary>
        /// Deletes row at specifies index.
        /// </summary>
        /// <param name="i">One-based index at which to delete.</param>
        public void DeleteRow(int i)
        {
            if (i <= 0 || i > rowCount)
                throw new ArgumentException("Index must be positive and <= number of rows.");

            Values.RemoveAt(i - 1);
            rowCount--;
        }

        /// <summary>
        /// Deletes column at specifies index.
        /// </summary>
        /// <param name="j">One-based index at which to delete.</param>
        public void DeleteColumn(int j)
        {
            if (j <= 0 || j > columnCount)
                throw new ArgumentException("Index must be positive and <= number of cols.");

            for (int i = 0; i < rowCount; i++)
            {
                ((ArrayList)Values[i]).RemoveAt(j - 1);
            }

            columnCount--;
        }

        /// <summary>
        /// Retrieves row vector at specfifed index and deletes it from matrix.
        /// </summary>
        /// <param name="i">One-based index at which to extract.</param>
        /// <returns>Row vector.</returns>
        public Matrix ExtractRow(int i)
        {
            Matrix buf = this.Row(i);
            this.DeleteRow(i);

            return buf;
        }

        /// <summary>
        /// Retrieves column vector at specfifed index and deletes it from matrix.
        /// </summary>
        /// <param name="j">One-based index at which to extract.</param>
        /// <returns>Row vector.</returns>
        public Matrix ExtractColumn(int j)
        {
            if (j <= 0 || j > columnCount)
                throw new ArgumentException("Index must be positive and <= number of cols.");

            Matrix buf = this.Column(j);
            this.DeleteColumn(j);

            return buf;
        }

        /// <summary>
        /// Inserts row at specified index.
        /// </summary>
        /// <param name="row">Vector to insert</param>
        /// <param name="i">One-based index at which to insert</param>
        public void InsertRow(Matrix row, int i)
        {
            int size = row.VectorLength();

            if (size == 0)
                throw new InvalidOperationException("Row must be a vector of length > 0.");

            if (i <= 0)
                throw new ArgumentException("Row index must be positive.");


            if (i > rowCount)
                this[i, size] = Complex.Zero;

            else if (size > columnCount)
            {
                this[i, size] = Complex.Zero;
                rowCount++;
            }
            else
                rowCount++;



            Values.Insert(--i, new ArrayList(size));
            //Debug.WriteLine(Values.Count.ToString());

            for (int k = 1; k <= size; k++)
            {
                ((ArrayList)Values[i]).Add(row[k]);
            }

            // fill w/ zeros if vector row is too short
            for (int k = size; k < columnCount; k++)
            {
                ((ArrayList)Values[i]).Add(Complex.Zero);
            }
        }

        /// <summary>
        /// Inserts a sub matrix M at row i and column j.
        /// </summary>
        /// <param name="i">One-based row number to insert.</param>
        /// <param name="j">One-based column number to insert.</param>
        /// <param name="M">Sub matrix to insert.</param>
        public void Insert(int i, int j, Matrix M)
        {            
            for (int m = 1; m <= M.rowCount; m++)
                for (int n = 1; n <= M.columnCount; n++)
                    this[i + m - 1, j + n - 1] = M[m, n];
        }

        /// <summary>
        /// Inserts column at specified index.
        /// </summary>
        /// <param name="col">Vector to insert</param>
        /// <param name="j">One-based index at which to insert</param>
        public void InsertColumn(Matrix col, int j)
        {
            int size = col.VectorLength();

            if (size == 0)
                throw new InvalidOperationException("Row must be a vector of length > 0.");

            if (j <= 0)
                throw new ArgumentException("Row index must be positive.");


            if (j > columnCount)
            {
                this[size, j] = Complex.Zero;
            }
            else
                columnCount++;

            if (size > rowCount)
            {
                this[size, j] = Complex.Zero;
            }

            j--;

            for (int k = 0; k < size; k++)
            {
                ((ArrayList)Values[k]).Insert(j, col[k + 1]);
            }

            // fill w/ zeros if vector col too short
            for (int k = size; k < rowCount; k++)
            {
                ((ArrayList)Values[k]).Insert(j, 0);
            }

        }

        /// <summary>
        /// Inverts square matrix as long as det != 0.
        /// </summary>
        /// <returns>Inverse of matrix.</returns>
        public Matrix Inverse()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot invert non-square matrix.");

            Complex det = this.Determinant();

            if (det == Complex.Zero)
                throw new InvalidOperationException("Cannot invert (nearly) singular matrix.");

            int n = this.columnCount;

            if (n == 1) return new Matrix(1 / det);

            if (this.IsReal() && this.IsOrthogonal())
                return this.Transpose();
            else if (this.IsUnitary())
                return this.ConjTranspose();

            if (this.IsDiagonal())
            {
                Matrix d = this.DiagVector();

                for (int i = 1; i <= n; i++)
                    d[i] = 1 / d[i];

                return Diag(d);
            }

            Complex[,] buf = new Complex[n, n];

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    buf[i, j] = Math.Pow(-1, i + j) * this.Minor(j + 1, i + 1).Determinant();
                }
            }

            return (new Matrix(buf) / det);
        }

        /// <summary>
        /// Alternative matrix inversion using Leverrier's formula
        /// </summary>
        /// <returns>Inverse of matrix.</returns>
        public Matrix InverseLeverrier()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot invert non-square matrix.");
            //else if (this.Determinant() == 0)
            //    throw new InvalidOperationException("Cannot invert (nearly) singular matrix.");

            int n = this.rowCount;
            Matrix Id = Identity(n);
            Matrix B = Id;
            Complex alpha;

            for (int k = 1; k < n; k++)
            {                
                Matrix buf = (this * B); // DEBUG                
                Complex buf2 = buf.Trace(); // DEBUG
                alpha = ((double)1 / k) * buf.Trace();
                B = alpha * Id - buf;
            }

            Matrix buf3 = (this * B); // DEBUG                
            Complex buf4 = buf3.Trace(); // DEBUG
            alpha = (this * B).Trace() / n;
            if (alpha != Complex.Zero)
                return B / alpha;
            else
                throw new InvalidOperationException("WARNING: Matrix nearly singular or badly scaled.");
        }

        /// <summary>
        /// Calcs the matrix that results in the clearing of a
        /// specified row and a specified column
        /// </summary>        
        /// <param name="row"></param>
        /// <param name="col"></param>
        /// <returns></returns>
        public Matrix Minor(int i, int j)
        {
            // THIS IS THE LOW-LEVEL SOLUTION ~ O(n^2)
            //Complex[,] buf = new Complex[RowCount - 1, ColumnCount - 1];
            //int r = 0;
            //int c = 0;

            //for (int i = 1; i <= RowCount; i++)
            //{
            //    if (i != row)
            //    {
            //        for (int j = 1; j <= ColumnCount; j++)
            //        {
            //            if (j != col)
            //            {
            //                buf[r, c] = this[i, j];
            //                c++;
            //            }
            //        }

            //        c = 0;
            //        r++;
            //    }
            //}

            //return new Matrix(buf);  

            // THIS IS THE HIGH-LEVEL SOLUTION ~ O(n)

            Matrix A = this.Clone();

            A.DeleteRow(i);
            A.DeleteColumn(j);

            return A;
        }

        /// <summary>
        /// Provides a shallow copy of this matrix in O(m).
        /// </summary>
        /// <returns></returns>
        public Matrix Clone()
        {
            Matrix A = new Matrix();
            A.rowCount = rowCount;
            A.columnCount = columnCount;

            for (int i = 0; i < rowCount; i++)
                A.Values.Add(((ArrayList)this.Values[i]).Clone());

            return A;
        }

        /// <summary>
        /// Extracts main diagonal vector of the matrix as a column vector.
        /// </summary>
        /// <returns></returns>
        public Matrix DiagVector()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot get diagonal of non-square matrix.");

            Matrix v = new Matrix(this.columnCount, 1);

            for (int i = 1; i <= this.columnCount; i++)
            {
                v[i] = this[i, i];
            }

            return v;
        }

        /// <summary>
        /// Retrieves column with one-based index j.
        /// </summary>
        /// <param name="j"></param>
        /// <returns>j-th column...</returns>
        public Matrix Column(int j)
        {
            Matrix buf = new Matrix(this.rowCount, 1);

            for (int i = 1; i <= this.rowCount; i++)
            {
                buf[i] = this[i, j];
            }

            return buf;
        }

        /// <summary>
        /// Retrieves row with one-based index i.
        /// </summary>
        /// <param name="i"></param>
        /// <returns>i-th row...</returns>
        public Matrix Row(int i)
        {
            if (i <= 0 || i > rowCount)
                throw new ArgumentException("Index exceed matrix dimension.");

            //return (new Matrix((Complex[])((ArrayList)Values[i - 1]).ToArray(typeof(Complex)))).Transpose();

            Matrix buf = new Matrix(columnCount, 1);

            for (int j = 1; j <= this.columnCount; j++)
            {
                buf[j] = this[i, j];
            }

            return buf;
        }

        /// <summary>
        /// Swaps each matrix entry A[i, j] with A[j, i].
        /// </summary>
        /// <returns>A transposed matrix.</returns>
        public Matrix Transpose()
        {
            Matrix M = new Matrix(columnCount, rowCount);

            for (int i = 1; i <= columnCount; i++)
            {
                for (int j = 1; j <= rowCount; j++)
                {
                    M[i, j] = this[j, i];
                }
            }

            return M;
        }

        /// <summary>
        /// Replaces each matrix entry z = x + iy with x - iy.
        /// </summary>
        /// <returns>Conjugated matrix.</returns>
        public Matrix Conjugate()
        {
            Matrix M = new Matrix(rowCount, columnCount);

            for (int i = 1; i <= rowCount; i++)
            {
                for (int j = 1; j <= columnCount; j++)
                {
                    M[i, j] = Complex.Conj(this[i, j]);
                }
            }

            return M;
        }

        /// <summary>
        /// Conjuagtes and transposes a matrix.
        /// </summary>
        /// <returns></returns>
        public Matrix ConjTranspose()
        {
            return this.Transpose().Conjugate();
        }

        /// <summary>
        /// Performs LU-decomposition of this instance and saves L and U
        /// within, where the diagonal elements belong to U
        /// (the ones of L are ones...)
        /// </summary>
        public void LU()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot perform LU-decomposition of non-square matrix.");

            int n = this.columnCount;

            for (int j = 1; j <= n; j++)
            {
                if (this[j, j] == 0)
                    throw new DivideByZeroException("Warning: Matrix badly scaled or close to singular. Try LUSafe() instead. Check if det != 0.");

                for (int k = 1; k < j; k++)
                {
                    for (int i = k + 1; i <= n; i++)
                    {
                        this[i, j] = this[i, j] - this[i, k] * this[k, j];
                    }
                }

                for (int i = j + 1; i <= n; i++)
                {
                    this[i, j] = this[i, j] / this[j, j];
                }
            }
        }

        /// <summary>
        /// Performs safe LU-decomposition of this instance with column pivoting 
        /// and saves L and U
        /// within, where the diagonal elements belong to U
        /// (the ones of L are ones...)
        /// </summary>
        /// <returns>Permutation matrix P with P*this = L*U</returns>
        /// <remarks>This needs additional time O(n^2).</remarks>
        public Matrix LUSafe()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot perform LU-decomposition of non-square matrix.");

            int n = this.columnCount;

            Matrix P = Identity(n); // permutation matrix
            int m;

            for (int j = 1; j <= n; j++)
            {
                //--> this test means probably deceleration
                //if (j < n && this.Extract(j + 1, n, j, j) == Zeros(n - j, 1))
                //    continue;

                #region Column pivoting

                // find index m with |this[m,j]| >= |this[i,j]| for all i in {j,...,n}
                if (j < n)
                {

                    m = j;

                    for (int i = j + 1; i <= n; i++)
                        if (Complex.Abs(this[i, j]) > Complex.Abs(this[m, j]))
                            m = i;

                    if (m > j) // <=> j2 != j
                    {
                        P.SwapRows(j, m);
                        this.SwapRows(j, m);
                    }

                    if (this[j, j] == 0)
                        throw new DivideByZeroException("Warning: Matrix close to singular.");
                }

                #endregion              

                for (int k = 1; k < j; k++)
                {
                    for (int i = k + 1; i <= n; i++)
                    {
                        this[i, j] = this[i, j] - this[i, k] * this[k, j];
                    }
                }

                for (int i = j + 1; i <= n; i++)
                {
                    this[i, j] = this[i, j] / this[j, j];
                }
            }

            return P;
        }


        /// <summary>
        /// Performs Cholesky decomposition of square, symmetric and positive definite
        /// matrix A = LL', where L is a lower triangular matrix. L is saved in the
        /// lower triangular part of A.</summary>
        /// <remarks>
        /// The diagonal elements can be retrieved
        /// by a_{11} = h_{11}^2, a_{ii} = h_{ii}^2 + \sum_{k=1}^{i-1}h_{ik}^2 (i = 2..n).
        /// Use CholeskyUndo() for convenience.
        /// WARNING: Cholesky decomposition only works for symmetric positive definite matrices!
        /// </remarks>        
        public void Cholesky()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot perform Cholesky decomposition of non-square matrix.");
           
            if (!this.IsSPD())
                throw new InvalidOperationException("Cannot perform Cholesky decomposition of matrix not being symmetric positive definite.");

            int n = rowCount;

            for (int k = 1; k < n; k++)
            {
                this[k, k] = Complex.Sqrt(this[k, k]);

                for (int i = 1; i <= n - k; i++)
                    this[k + i, k] = this[k + i, k] / this[k, k];

                for (int j = k + 1; j <= n; j++)
                    for (int i = 0; i <= n - j; i++)
                        this[j + i, j] = this[j + i, j] - this[j + i, k] * this[j, k];
            }

            this[n, n] = Complex.Sqrt(this[n, n]);

        }

        /// <summary>
        /// Since the cholesky decomposition is saved within the symmetric matrix to be
        /// decomposited, it can be undone to restore the initial matrix.
        /// </summary>
        public void CholeskyUndo()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot undo cholesky decomposition on non-square matrix.");

            this[1, 1] = Sqr(this[1, 1]);

            Complex buf;

            for (int i = 2; i <= rowCount; i++)
            {
                buf = Complex.Zero;

                for (int k = 1; k <= i - 1; k++)
                    buf += Sqr(this[i, k]);

                this[i, i] = Sqr(this[i, i]) + buf;
            }

            this.SymmetrizeDown();
        }

        /// <summary>
        /// Performs forward insertion for regular lower triangular matrix
        /// and right side b, such that the solution is saved right within b.
        /// The matrix is not changed.
        /// </summary>
        /// <param name="b">Vector of height n, if matrix is n by n.</param>
        public void ForwardInsertion(Matrix b)
        {
            if (!this.IsLowerTriangular())
                throw new InvalidOperationException("Cannot perform forward insertion for matrix not being lower triangular.");

            if (/*this.Determinant*/this.DiagProd() == 0)
                throw new InvalidOperationException("Warning: Matrix is nearly singular.");

            int n = rowCount;

            if (b.VectorLength() != n)
                throw new ArgumentException("Parameter must vector of the same height as matrix.");

            for (int j = 1; j <= n - 1; j++)
            {
                b[j] /= this[j, j];

                for (int i = 1; i <= n - j; i++)
                    b[j + i] -= b[j] * this[j + i, j];
            }

            b[n] /= this[n, n];
        }

        /// <summary>
        /// Performs backward insertion for regular upper triangular matrix
        /// and right side b, such that the solution is saved right within b.
        /// The matrix is not changed.
        /// </summary>
        /// <param name="b">Vector of height n, if matrix is n by n.</param>
        public void BackwardInsertion(Matrix b)
        {
            if (!this.IsUpperTriangular())
                throw new InvalidOperationException("Cannot perform backward insertion for matrix not being upper triangular.");

            if (/*this.Determinant*/this.DiagProd() == 0)
                throw new InvalidOperationException("Warning: Matrix is nearly singular.");

            int n = rowCount;

            if (b.VectorLength() != n)
                throw new ArgumentException("Parameter must vector of the same height as matrix.");

            for (int j = n; j >= 2; j--)
            {
                b[j] /= this[j, j];

                for (int i = 1; i <= j - 1; i++)
                    b[i] -= b[j] * this[i, j];
            }

            b[1] /= this[1, 1];
        }
       
        /// <summary>
        /// Makes square matrix symmetric by copying the upper half to the lower half.
        /// </summary>
        public void SymmetrizeDown()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot symmetrize non-square matrix.");

            for (int j = 1; j <= columnCount; j++)
                for (int i = j + 1; i <= columnCount; i++)
                    this[i, j] = this[j, i];
        }

        /// <summary>
        /// Makes square matrix symmetric by copying the lower half to the upper half.
        /// </summary>
        public void SymmetrizeUp()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot symmetrize non-square matrix.");

            for (int i = 1; i <= rowCount; i++)
                for (int j = i + 1; j <= columnCount; j++)
                    this[i, j] = this[j, i];
        }

        /// <summary>
        /// Gram-Schmidtian orthogonalization of an m by n matrix A, such that
        /// {Q, R} is returned, where A = QR, Q is m by n and orthogonal, R is
        /// n by n and upper triangular matrix.
        /// </summary>
        /// <returns></returns>
        public Matrix[] QRGramSchmidt()
        {
            int m = rowCount;
            int n = columnCount;

            Matrix A = this.Clone();

            Matrix Q = new Matrix(m, n);
            Matrix R = new Matrix(n, n);

            // the first column of Q equals the first column of this matrix
            for (int i = 1; i <= m; i++)
                Q[i, 1] = A[i, 1];

            R[1, 1] = Complex.One;

            for (int k = 1; k <= n; k++)
            {
                R[k, k] = new Complex(A.Column(k).Norm());

                for (int i = 1; i <= m; i++)
                    Q[i, k] = A[i, k] / R[k, k];

                for (int j = k + 1; j <= n; j++)
                {
                    R[k, j] = Dot(Q.Column(k), A.Column(j));

                    for (int i = 1; i <= m; i++)
                        A[i, j] = A[i, j] - Q[i, k] * R[k, j];
                }
            }

            return new Matrix[] { Q, R };
        }

        /// <summary>
        /// Computes approximates of the eigenvalues of this matrix. WARNING: Computation
        /// uses basic QR iteration with Gram-Schmidtian orthogonalization. This implies that
        /// (1) only real matrices can be examined; (2) if the matrix has a multiple eigenvalue
        /// or complex eigenvalues, partial junk is returned. This is due to the eigenvalues having
        /// to be like |L1| > |L2| > ... > |Ln| for QR iteration to work properly.
        /// </summary>
        /// <returns></returns>
        public Matrix Eigenvalues()
        {
            return this.QRIterationBasic(40).DiagVector();

        }

        /// <summary>
        /// Computes eigenvector from eigenvalue.
        /// </summary>
        /// <param name="eigenvalue"></param>
        /// <returns></returns>
        public Matrix Eigenvector(Complex eigenvalue)
        {

            throw new NotImplementedException();
            
        }

        /// <summary>
        /// Solves equation this*x = b via conjugate gradient method.
        /// </summary>
        /// <param name="b"></param>
        /// <returns></returns>
        public Matrix SolveCG(Matrix b)
        {
            throw new NotImplementedException("Still buggy!");

            if (!this.IsSPD())
                throw new InvalidOperationException("CG method only works for spd matrices.");
            else if (!this.IsReal())
                throw new InvalidOperationException("CG method only works for real matrices.");

            int n = rowCount;
            int max_iterations = 150;
            double tolerance = 1e-6;

            Matrix x = Ones(n, 1); // x will contain the solution
            Matrix r = b - this * x; // residual approaches zero as x converges to the solution
            Matrix d = r; // dir = direction of descence
            double delta = r.Norm(); // delta denotes the current error
            delta *= delta; 
            tolerance *= tolerance;

            Matrix h = Zeros(n, 1);
            double alpha, gamma;
            double old_delta;

            if (delta <= tolerance)
                return x;
            else
            {
                for (int i = 0; i < max_iterations; i++)
                {
                    h = this * d;
                    gamma = Dot(h, d).Re;

                    if (Math.Abs(gamma) <= tolerance)
                        return x;

                    alpha = delta / gamma;

                    x += alpha * d; // compute new approximation of solution
                    r -= alpha * h; // compute new residual

                    old_delta = delta; // buffer delta

                    delta = r.Norm();
                    delta *= delta;

                    if (delta <= tolerance)
                        return x;

                    d = r + delta / old_delta * d; // compute new direction of descence
                }

                return x;
            }
        }

        /// <summary>
        /// Executes the QR iteration.
        /// </summary>
        /// <param name="max_iterations"></param>
        /// <returns></returns>
        public Matrix QRIterationBasic(int max_iterations)
        {
            if (!this.IsReal())
                throw new InvalidOperationException("Basic QR iteration is possible only for real matrices.");

            Matrix T = this.Clone();
            Matrix[] QR = new Matrix[2];

            for (int i = 0; i < max_iterations; i++)
            {
                QR = T.QRGramSchmidt();
                T = QR[1] * QR[0];
            }

            return T;
        }

        /// <summary>
        /// QR iteration using Hessenberg-Householder reduction.
        /// </summary>
        /// <param name="max_iterations"></param>
        /// <returns></returns>
        public Matrix QRIterationHessenberg(int max_iterations)
        {

            throw new NotImplementedException("Still buggy!");

            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot perform QR iteration of non-square matrix.");

            int n = this.RowCount;

            Matrix[] TQ = this.HessenbergHouseholder();
            Matrix T = TQ[0];

            for (int j = 1; j <= max_iterations; j++)
            {
                Matrix[] QRcs = T.QRGivens();
                T = QRcs[1];

                for (int k = 1; k <= n - 1; k++)
                {
                    T.Gacol(QRcs[2][k], QRcs[3][k], 1, k + 1, k, k + 1);
                }
            }

            return T;
        }

        /// <summary>
        /// QR factorization avec Givens rotations.
        /// </summary>
        /// <param name="H"></param>
        /// <returns></returns>
        public Matrix[] QRGivens()
        {

            throw new NotImplementedException("Still buggy!");

            Matrix H = this.Clone();
            int m = H.RowCount;
            int n = H.ColumnCount;

            Matrix c = Zeros(n - 1, 1);
            Matrix s = Zeros(n - 1, 1);
            Complex[] cs;

            for (int k = 1; k <= n - 1; k++)
            {
                cs = GivensCS(H[k, k], H[k + 1, k]);
                c[k] = cs[0];
                s[k] = cs[1];
                this.Garow(c[k], s[k], 1, k + 1, k, k + 1);
            }

            return new Matrix[] { GivProd(c, s, n), H, c, s };
        }

        /// <summary>
        /// Givens product. Internal use for QRGivens.
        /// </summary>
        /// <param name="c"></param>
        /// <param name="s"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private Matrix GivProd(Matrix c, Matrix s, int n)
        {
            int n1 = n - 1;
            int n2 = n - 2;

            Matrix Q = Eye(n);
            Q[n1, n1] = c[n1];
            Q[n, n] = c[n1];
            Q[n1, n] = s[n1];
            Q[n, n1] = -s[n1];

            for (int k = n2; k >= 1; k--)
            {
                int k1 = k + 1;
                Q[k, k] = c[k];
                Q[k1, k] = -s[k];
                Matrix q = Q.Extract(k1, k1, k1, n);
                Q.Insert(k, k1, s[k] * q);
                Q.Insert(k1, k1, c[k] * q);
            }

            return Q;
        }

        /// <summary>
        /// Product G(i,k,theta)'*this. Internal use for QRGivens.
        /// </summary>
        /// <param name="c"></param>
        /// <param name="s"></param>
        /// <param name="i"></param>
        /// <param name="k"></param>
        /// <param name="j1"></param>
        /// <param name="j2"></param>
        /// <returns></returns>
        private void Garow(Complex c, Complex s, int i, int k, int j1, int j2)
        {
            for (int j = j1; j <= j2; j++)
            {
                Complex t1 = this[i, j];
                Complex t2 = this[k, j];
                this[i, j] = c * t1 - s * t2;
                this[k, j] = s * t1 + c * t2;
            }
        }

        /// <summary>
        /// Product M*G(i,k,theta). Internal use for QRGivens.
        /// </summary>
        /// <param name="c"></param>
        /// <param name="s"></param>
        /// <param name="j1"></param>
        /// <param name="j2"></param>
        /// <param name="i"></param>
        /// <param name="k"></param>
        public void Gacol(Complex c, Complex s, int j1, int j2, int i, int k)
        {
            for (int j = j1; j <= j2; j++)
            {
                Complex t1 = this[j, i];
                Complex t2 = this[j, k];

                this[j, i] = c * t1 - s * t2;
                this[j, k] = s * t1 + c * t2;
            }
        }

        /// <summary>
        /// Computes Givesn sine and cosine.
        /// </summary>
        /// <param name="xi"></param>
        /// <param name="xk"></param>
        /// <returns></returns>
        private Complex[] GivensCS(Complex xi, Complex xk)
        {
            Complex c = Complex.Zero;
            Complex s = Complex.Zero;

            if (xk == 0)
            {
                c = Complex.One;
            }
            else
            {
                if (Complex.Abs(xk) > Complex.Abs(xi))
                {
                    Complex t = -xi / xk;
                    s = 1 / (Complex.Sqrt(1 + t * t));
                    c = s * t;
                }
                else
                {
                    Complex t = -xk / xi;
                    c = 1 / (Complex.Sqrt(1 + t * t));
                    s = c * t;
                }
            }

            return new Complex[] { c, s };
        }

        #endregion

        #region Numbers

        /// <summary>
        /// Calcs determinant of square matrix
        /// </summary>
        /// <returns></returns>
        public Complex Determinant()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot calc determinant of non-square matrix.");
            else if (this.columnCount == 1)
                return this[1, 1];
            else if (this.IsTrapeze()) // is square, therefore triangular
            {
                return this.DiagProd();
            }
            else
            {
                // perform LU-decomposition & return product of diagonal elements of U
                Matrix X = this.Clone();

                // for speed concerns, use this
                //X.LU();
                //return X.DiagProd();

                // this is slower and needs more memory... .
                Matrix P = X.LUSafe();
                return (double)P.Signum() * X.DiagProd();
            }
        }

        public double ColumnSumNorm()
        {
            return this.TaxiNorm();
        }

        public double RowSumNorm()
        {
            return this.MaxNorm();
        }

        /// <summary>
        /// Computes the permanent of the current instance. WARNING: This algorithm has exponential runtime.
        /// Don't use for any but very small instances.
        /// </summary>
        /// <returns></returns>
        public Complex Permanent()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot compute permanent of non-square matrix.");

            if (this.HasZeroRowOrColumn())
                return Complex.Zero;

            if (this == Ones(rowCount))
                return new Complex(Factorial(rowCount));

            if (this.IsPermutation())
                return Complex.One;            

            Complex buf = Complex.Zero;

            int minRow = this.GetMinRow();
            int minColumn = this.GetMinColumn();

            if (this.AbsRowSum(minRow) < this.AbsColumnSum(minColumn))
            {
                for (int j = 1; j <= columnCount; j++)
                    if (this[minRow, j] != 0)
                        buf += this[minRow, j] * this.Minor(minRow, j).Permanent();
            }
            else
            {
                for (int i = 1; i <= rowCount; i++)
                    if (this[i, minColumn] != 0)
                        buf += this[i, minColumn] * this.Minor(i, minColumn).Permanent();
            }

            return buf;
        }

        /// <summary>
        /// Finds index of row with minimal AbsRowSum.
        /// </summary>
        /// <returns></returns>
        public int GetMinRow()
        {
            double buf = this.AbsRowSum(1);
            int index = 1;

            double buf2;

            for (int i = 2; i <= rowCount; i++)
            {
                buf2 = this.AbsRowSum(i);
                if (buf2 < buf)
                {
                    buf = buf2;
                    index = i;
                }
            }

            return index;
        }

        /// <summary>
        /// Finds index of column with minimal AbsColumnSum.
        /// </summary>
        /// <returns></returns>
        public int GetMinColumn()
        {
            double buf = this.AbsColumnSum(1);
            int index = 1;

            double buf2;

            for (int j = 2; j <= columnCount; j++)
            {
                buf2 = this.AbsColumnSum(j);
                if (buf2 < buf)
                {
                    buf = buf2;
                    index = j;
                }
            }

            return index;
        }



        private double Factorial(int x)
        {
            double buf = 1;
            for (int i = 2; i <= x; i++)
                buf *= i;

            return buf;
        }

        /// <summary>
        /// Computes signum of a permutation matrix, which is 1 for an even
        /// number of swaps and -1 for an odd number of swaps. WARNING: 
        /// if *this is not a permutation matrix (e.i. a permutation of Id),
        /// garbage is returned.
        /// </summary>
        /// <returns></returns>
        public double Signum()
        {
            double buf = 1;

            int n = rowCount;
            double fi, fj;

            for (int i = 1; i < n; i++)
            {
                for (fi = 1; fi < n && this[i, (int)fi] != Complex.One; fi++) ;

                for (int j = i + 1; j <= n; j++)
                {
                    for (fj = 1; fj <= n && this[j, (int)fj] != Complex.One; fj++) ;

                    buf *= (fi - fj) / (i - j);
                }
            }

            return buf;
        }     

        /// <summary>
        ///  Calcs condition number with respect to inversion
        /// by using |A|*|inv(A)| and 1-Norm.
        /// </summary>
        /// <returns></returns>
        public double Condition()
        {
            return this.TaxiNorm() * this.Inverse().TaxiNorm();
        }

        /// <summary>
        ///  Calcs condition number with respect to inversion
        /// by using |A|*|inv(A)| and p norm.
        /// </summary>
        /// <param name="p">Specifies the norm to be used. Can be one or positive infinity.</param>
        /// <returns></returns>
        public double Condition(int p)
        {
            return this.PNorm(p) * this.Inverse().PNorm(p);
        }

        /// <summary>
        ///  Calcs condition number with respect to inversion
        /// by using |A|*|inv(A)| and frobenius norm.
        /// </summary>
         /// <returns></returns>
        public double ConditionFro()
        {
            return this.FrobeniusNorm() * this.Inverse().FrobeniusNorm();
        }

        /// <summary>
        /// Calcs p-norm of given matrix: p-th root of the sum
        /// of the p-th powers of the absolute values of all matrix entries.
        /// </summary>
        /// <param name="p">Which norm to compute; can be positive infinity.</param>
        /// <returns></returns>
        /// <remarks>If p not in {i, +infty}, *this must be a vector.</remarks>
        public double PNorm(double p)
        {
            if (p <= 0)
                throw new ArgumentException("Argument must be greater than zero.");

            if (p == 1) return TaxiNorm();
            else if (p == double.PositiveInfinity) return MaxNorm();

            int dim = this.VectorLength();
            if (dim == 0)
                throw new InvalidOperationException("Cannot calc p-norm of matrix.");

            double buf = 0;

            for (int i = 1; i <= dim; i++)
                buf += Math.Pow(Complex.Abs(this[i]), p);

            return Math.Pow(buf, 1 / p);
        }

        /// <summary>
        /// 2-Norm for vectors. If *this is a matrix, you might want to choose
        /// FrobeniusNorm().
        /// </summary>
        /// <returns></returns>
        public double Norm()
        {
            return PNorm(2);
        }

        /// <summary>
        /// Frobenius norm of a square matrix. If *this is a vector, this method
        /// is equivalent to Norm() and PNorm(2).
        /// </summary>
        /// <returns></returns>
        public double FrobeniusNorm()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot compute frobenius norm of non-square matrix.");

            int n = this.columnCount;
            double buf = 0;

            for (int i = 1; i <= n; i++)
                for (int j = 1; j <= n; j++)
                    buf += (this[i, j] * Complex.Conj(this[i, j])).Re;

            return Math.Sqrt(buf);
        }
  
        /// <summary>
        /// Also known as column-sum norm.
        /// </summary>
        /// <returns>Maximal AbsColumnSum</returns>
        public double TaxiNorm()
        {
            double buf = 0;

            int dim = this.VectorLength();

            if (dim != 0) // vector case
            {
                for (int i = 1; i <= dim; i++)
                {
                    buf += Complex.Abs(this[i]);
                }
            }
            else // general case
            {
                double buf2 = 0;

                for (int j = 1; j <= columnCount; j++)
                {
                    buf2 = AbsColumnSum(j);
                    if (buf2 > buf)
                        buf = buf2;
                }
            }

            return buf;
        }

        /// <summary>
        /// Also known as row-sum norm.
        /// </summary>
        /// <returns>Maximal AbsRowSum</returns>
        public double MaxNorm()
        {
            double buf = 0;
            double buf2 = 0;

            int dim = this.VectorLength();

            if (dim != 0) // vector case
            {
                for (int i = 1; i <= dim; i++)
                {
                    buf2 = Complex.Abs(this[i]);
                    if (buf2 > buf)
                        buf = buf2;
                }
            }
            else // general case
            {
                for (int i = 1; i <= rowCount; i++)
                {
                    buf2 = AbsRowSum(i);
                    if (buf2 > buf)
                        buf = buf2;
                }
            }

            return buf;
        }

        /// <summary>
        /// Calcs sum of the elements of a certain col.
        /// </summary>
        /// <param name="i">One-based index of the col to consider.</param>
        /// <returns></returns>
        public Complex ColumnSum(int j)
        {
            if (j <= 0 || j > columnCount)
                throw new ArgumentException("Index out of range.");

            Complex buf = Complex.Zero;

            j--;

            for (int i = 0; i < rowCount; i++)
            {
                buf += (Complex)(((ArrayList)Values[i])[j]);
            }

            return buf;
        }

        /// <summary>
        /// Calcs sum of the absolute values of the elements of a certain col.
        /// </summary>
        /// <param name="i">One-based index of the col to consider.</param>
        /// <returns></returns>
        public double AbsColumnSum(int j)
        {
            if (j <= 0 || j > columnCount)
                throw new ArgumentException("Index out of range.");

            double buf = 0;
           
            for (int i = 1; i <= rowCount; i++)
            {
                buf += Complex.Abs(this[i,j]);
            }

            return buf;
        }

        /// <summary>
        /// Calcs sum of the elements of a certain row.
        /// </summary>
        /// <param name="i">One-based index of the row to consider.</param>
        /// <returns></returns>
        public Complex RowSum(int i)
        {
            if (i <= 0 || i > rowCount)
                throw new ArgumentException("Index out of range.");

            Complex buf = Complex.Zero;

            ArrayList row = (ArrayList)Values[i - 1];

            for (int j = 0; j < columnCount; j++)
            {
                buf += (Complex)(row[j]);
            }

            return buf;
        }

        /// <summary>
        /// Calcs sum of the absolute values of the elements of a certain row.
        /// </summary>
        /// <param name="i">One-based index of the row to consider.</param>
        /// <returns></returns>
        public double AbsRowSum(int i)
        {
            if (i <= 0 || i > rowCount)
                throw new ArgumentException("Index out of range.");

            double buf = 0;            

            for (int j = 1; j <= columnCount; j++)
            {
                buf += Complex.Abs(this[i, j]);
            }

            return buf;
        }

        /// <summary>
        /// Computes product of main diagonal entries.
        /// </summary>
        /// <returns>Product of diagonal elements</returns>
        public Complex DiagProd()
        {
            Complex buf = Complex.One;
            int dim = Math.Min(this.rowCount, this.columnCount);

            for (int i = 1; i <= dim; i++)
            {
                buf *= this[i, i];
            }

            return buf;
        }

        /// <summary>
        /// Calcs trace of the matrix.
        /// </summary>
        /// <returns>Sum of diagonal elements.</returns>
        public Complex Trace()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot calc trace of non-square matrix.");

            Complex buf = Complex.Zero;

            for (int i = 1; i <= this.rowCount; i++)
            {
                buf += this[i, i];
            }

            return buf;
        }

        Complex Sqr(Complex x)
        {
            return x * x;
        }

        #endregion

        #region Checks

        /// <summary>
        /// Matrix is normal, iff A*A^H = A^H*A, where A is the conjugated transposed of A.
        /// </summary>
        /// <returns></returns>
        public bool IsNormal()
        {
            return (this * this.ConjTranspose() == this.ConjTranspose() * this);
        }

        /// <summary>
        /// Matrix is unitary, iff A^H*A = id, where A^H is the conjugated transpose of A.
        /// </summary>
        /// <returns>True iff matrix is unitary.</returns>
        public bool IsUnitary()
        {
            if(!this.IsSquare())
                return false;

            return (this.ConjTranspose() * this == Identity(rowCount));
        }

        /// <summary>
        /// Matrix A is Hermitian iff A^H = A, where A^H is the conjugated transposed of A.
        /// </summary>
        /// <returns>True iff matrix is Hermitian.</returns>
        public bool IsHermitian()
        {
            if (!this.IsSquare())
                return false;

            return this.ConjTranspose() == this;
        }

        /// <summary>
        /// Checks if matrix consists only of real entries.
        /// </summary>
        /// <returns>True iff all entries are real.</returns>
        public bool IsReal()
        {
            for (int i = 1; i <= rowCount; i++)            
                for (int j = 1; j <= columnCount; j++)                
                    if (!this[i, j].IsReal()) return false;

            return true;
        }
        
        /// <summary>
        /// Checks for symmetric positive definiteness.
        /// </summary>
        /// <returns>True iff matrix is symmetrix positive definite.</returns>
        public bool IsSymmetricPositiveDefinite()
        {
            return (this.IsSymmetric() && this.Definiteness() == DefinitenessType.PositiveDefinite);
        }

        /// <summary>
        /// Checks for symmetric positive definiteness.
        /// </summary>
        /// <returns>True iff matrix is symmetrix positive definite.</returns>
        public bool IsSPD()
        {
            return this.IsSymmetricPositiveDefinite();
        }

        /// <summary>
        /// Finds out the type of definiteness of a symmetric square matrix.
        /// </summary>
        /// <returns></returns>
        public DefinitenessType Definiteness()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Definiteness undefined for non-square matrices.");
            else if (this == Zeros(this.rowCount, this.columnCount))
                return DefinitenessType.Indefinite;
            else if (!this.IsSymmetric())
                throw new InvalidOperationException("This test works only for symmetric matrices.");
            else if (!this.IsReal())
                throw new InvalidOperationException("This test only works for real matrices.");

            // step 1: construct orthogonal basis for A
            // using Gram-Schmidt orthogonalization
            int n = this.rowCount;
           
            Matrix[] y = new Matrix[n + 1];
            for (int i = 0; i <= n; i++)
                y[i] = Matrix.Zeros(n, 1);

            y[1] = this.Column(1);

            Matrix xk; // to buffer this.Column(k)
            Matrix buf;

            // Gram-Schmidt:
            for (int k = 2; k <= n; k++)
            {
                xk = this.Column(k);

                buf = Zeros(n, 1);

                for (int i = 1; i < k; i++)
                {
                    buf += y[i] * Dot(xk, this * y[i]) / Dot(y[i], this * y[i]);
                }

                y[k] = xk - buf;
            }

          

            // step 2: test for definiteness; 
            // e.g. A pos. def. <=> A > 0 <=> y[i]'Ay[i] > 0 for all i (same for neg. def., ...)

            bool strict = true; // pos. def || neg. def.
            Complex res;

            for (int i = 1; i < n; i++)
            {
                res = Dot(y[i], this * y[i]) * Dot(y[i + 1], this * y[i + 1]);

                if (res == 0) strict = false;
                else if (res.Re < 0) return DefinitenessType.Indefinite;
            }

            if (Dot(y[1], this * y[1]).Re >= 0)
            {
                if (strict) return DefinitenessType.PositiveDefinite;
                else return DefinitenessType.PositiveSemidefinite;
            }
            else
            {
                if (strict) return DefinitenessType.NegativeDefinite;
                else return DefinitenessType.NegativeSemidefinite;
            }
        }

        /// <summary>
        /// Checks if matrix has a row or column consisting of zeros.
        /// </summary>
        /// <returns>True iff so.</returns>
        public bool HasZeroRowOrColumn()
        {
            for (int i = 1; i <= rowCount; i++)
                if (this.AbsRowSum(i) == 0) return true;

            for (int i = 1; i <= columnCount; i++)
                if (this.AbsColumnSum(i) == 0) return true;

            return false;
        }

        /// <summary>
        /// Checks if matrix consists only of zeros and ones.
        /// </summary>
        /// <returns></returns>
        public bool IsZeroOneMatrix()
        {
            for (int i = 1; i <= rowCount; i++)
                for (int j = 1; j <= columnCount; j++)
                    if (this[i, j] != Complex.Zero && this[i, j] != Complex.One)
                        return false;

            return true;
        }

        /// <summary>
        /// Checks if matrix is permutation of the identity matrix.
        /// </summary>
        /// <returns>True iff matrix is permutation matrix.</returns>
        public bool IsPermutation()
        {
            return (!this.IsSquare() && this.IsZeroOneMatrix() && this.IsInvolutary());
                
        }

        /// <summary>
        /// Checks if matrix is diagonal matrix.
        /// </summary>
        /// <returns>True iff matrix is diagonal.</returns>
        public bool IsDiagonal()
        {
            return (this.Clone() - Diag(this.DiagVector()) == Zeros(rowCount, columnCount));
        }

        /// <summary>
        /// Checks if matrix is n by one or one by n.
        /// </summary>
        /// <returns>Length, if vector; zero else.</returns>
        public int VectorLength()
        {
            if (columnCount > 1 && rowCount > 1)
                return 0;
            else return Math.Max(columnCount, rowCount);
        }

        /// <summary>
        /// Checks if number of rows equals number of columns.
        /// </summary>
        /// <returns>True iff matrix is n by n.</returns>
        public bool IsSquare()
        {
            return (this.columnCount == this.rowCount);
        }

        /// <summary>
        /// Checks if matrix is involutary, e.i. if A*A = id.
        /// </summary>
        /// <returns>True iff matrix is involutary.</returns>
        public bool IsInvolutary()
        {
            return (this * this == Identity(rowCount));
        }

        /// <summary>
        /// Checks if A[i, j] == A[j, i].
        /// </summary>
        /// <returns>True iff matrix is symmetric.</returns>
        public bool IsSymmetric()
        {
            for (int i = 1; i <= this.rowCount; i++)
            {
                for (int j = 1; j <= this.columnCount; j++)
                {
                    if (this[i, j] != this[j, i])
                        return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Checks for orthogonality by testing if A*A' == id.
        /// </summary>
        /// <returns>True iff matrix is orthogonal.</returns>
        public bool IsOrthogonal()
        {
            return (this.IsSquare() && this * this.Transpose() == Identity(this.rowCount));
        }

        /// <summary>
        /// Checks if matrix is lower or upper trapeze.
        /// </summary>
        /// <returns>True iff matrix is trapeze.</returns>
        public bool IsTrapeze()
        {
            return (this.IsUpperTrapeze() || this.IsLowerTrapeze());
        }

        /// <summary>
        /// Checks if matrix is trapeze and square.
        /// </summary>
        /// <returns>True iff matrix is triangular.</returns>
        public bool IsTriangular()
        {
            return (this.IsLowerTriangular() || this.IsUpperTriangular());
        }

        /// <summary>
        /// Checks if matrix is square and upper trapeze.
        /// </summary>
        /// <returns>True iff matrix is upper triangular.</returns>
        public bool IsUpperTriangular()
        {
            return (this.IsSquare() && this.IsUpperTrapeze());
        }

        /// <summary>
        /// Checks if matrix is square and lower trapeze.
        /// </summary>
        /// <returns>True iff matrix is lower triangular.</returns>
        public bool IsLowerTriangular()
        {
            return (this.IsSquare() && this.IsLowerTrapeze());
        }

        /// <summary>
        /// Checks if A[i, j] == 0 for i < j.
        /// </summary>
        /// <returns>True iff matrix is upper trapeze.</returns>
        public bool IsUpperTrapeze()
        {
            for (int j = 1; j <= columnCount; j++)
                for (int i = j + 1; i <= rowCount; i++)
                    if (this[i, j] != 0) return false;

            return true;
        }

        /// <summary>
        /// Checks if A[i, j] == 0 for i > j.
        /// </summary>
        /// <returns>True iff matrix is lower trapeze.</returns>
        public bool IsLowerTrapeze()
        {
            for (int i = 1; i <= rowCount; i++)
                for (int j = i + 1; j <= columnCount; j++)
                    if (this[i, j] != 0) return false;

            return true;
        }

        #endregion       

        #endregion                

        #region Overrides & Operators

        public override string ToString()
        {
            string s = "";
            Complex buf;

            for (int i = 1; i <= rowCount; i++)
            {
                for (int j = 1; j <= columnCount; j++)
                {
                    buf = this[i, j];
                    if (buf.Re == double.PositiveInfinity || buf.Re == double.NegativeInfinity
                        || buf.Im == double.PositiveInfinity || buf.Im == double.NegativeInfinity)
                        s += "oo";                    
                    else if (buf.Re == double.NaN || buf.Im == double.NaN)
                        s += "?";
                    else
                        s += buf.ToString();

                    s += ";" + "\t";
                }

                s += "\\" + System.Environment.NewLine;
            }

            return s;
        }

        public string ToString(string format)
        {
            string s = "";
            Complex buf;

            for (int i = 1; i <= rowCount; i++)
            {
                for (int j = 1; j <= columnCount; j++)
                {
                    buf = this[i, j];
                    if (buf.Re == double.PositiveInfinity || buf.Re == double.NegativeInfinity
                        || buf.Im == double.PositiveInfinity || buf.Im == double.NegativeInfinity)
                        s += "oo";
                    else if (buf.Re == double.NaN || buf.Im == double.NaN)
                        s += "?";
                    else
                        s += buf.ToString(format);

                    s += ";" + "\t";
                }

                s += "\\" + System.Environment.NewLine;
            }

            return s;
        }

        public override bool Equals(object obj)
        {
            return obj.ToString() == this.ToString();
        }

        public override int GetHashCode()
        {
            return -1;
        }

        public static bool operator ==(Matrix A, Matrix B)
        {
           
            if (A.RowCount != B.RowCount || A.ColumnCount != B.ColumnCount)
                return false;

            for (int i = 1; i <= A.RowCount; i++)
            {
                for (int j = 1; j <= A.ColumnCount; j++)
                {
                    if (A[i, j] != B[i, j]) return false;
                }
            }

            return true;
        }

        public static bool operator !=(Matrix A, Matrix B)
        {
            return !(A == B);
        }

        public static Matrix operator +(Matrix A, Matrix B)
        {
          
            if (A.RowCount != B.RowCount || A.ColumnCount != B.ColumnCount)
                throw new ArgumentException("Matrices must be of the same dimension.");

            for (int i = 1; i <= A.RowCount; i++)
            {
                for (int j = 1; j <= A.ColumnCount; j++)
                {
                    A[i, j] += B[i, j];
                }
            }

            return A;
        }

        public static Matrix operator -(Matrix A, Matrix B)
        {
          
            if (A.RowCount != B.RowCount || A.ColumnCount != B.ColumnCount)
                throw new ArgumentException("Matrices must be of the same dimension.");

            for (int i = 1; i <= A.RowCount; i++)
            {
                for (int j = 1; j <= A.ColumnCount; j++)
                {
                    A[i, j] -= B[i, j];
                }
            }

            return A;
        }

        public static Matrix operator -(Matrix A)
        {
           
            for (int i = 1; i <= A.RowCount; i++)
            {
                for (int j = 1; j <= A.ColumnCount; j++)
                {
                    A[i, j] = -A[i, j];
                }
            }

            return A;
        }

        public static Matrix operator *(Matrix A, Matrix B)
        {
            
            if (A.ColumnCount != B.RowCount)
                throw new ArgumentException("Inner matrix dimensions must agree.");

            Matrix C = new Matrix(A.RowCount, B.ColumnCount);

            for (int i = 1; i <= A.RowCount; i++)
            {
                for (int j = 1; j <= B.ColumnCount; j++)
                {
                    C[i, j] = Dot(A.Row(i), B.Column(j));
                }
            }

            return C;

        }

        public static Matrix operator *(Matrix A, Complex x)
        {
           
            Matrix B = new Matrix(A.rowCount, A.columnCount);

            for (int i = 1; i <= A.RowCount; i++)
            {
                for (int j = 1; j <= A.ColumnCount; j++)
                {
                    B[i, j] = A[i, j] * x;
                }
            }

            return B;
        }

        public static Matrix operator *(Complex x, Matrix A)
        {
           
            Matrix B = new Matrix(A.RowCount, A.ColumnCount);

            for (int i = 1; i <= A.RowCount; i++)
            {
                for (int j = 1; j <= A.ColumnCount; j++)
                {
                    B[i, j] = A[i, j] * x;
                }
            }

            return B;
        }

        public static Matrix operator *(Matrix A, double x)
        {
            return (new Complex(x)) * A;
        }

        public static Matrix operator *(double x, Matrix A)
        {
            return (new Complex(x)) * A;
        }        

        public static Matrix operator /(Matrix A, Complex x)
        {
            return (1 / x) * A;
        }

        public static Matrix operator /(Matrix A, double x)
        {
            return (new Complex(1 / x)) * A;
        }

        public static Matrix operator ^(Matrix A, int k)
        {
            if (k < 0)
                if (A.IsSquare())
                    return A.InverseLeverrier() ^ (-k);
                else throw new InvalidOperationException("Cannot take non-square matrix to the power of zero.");
            else if (k == 0)
                if (A.IsSquare())
                    return Matrix.Identity(A.RowCount);
                else throw new InvalidOperationException("Cannot take non-square matrix to the power of zero.");
            else if (k == 1)
                if (A.IsSquare())
                    return A;
                else throw new InvalidOperationException("Cannot take non-square matrix to the power of one.");
            else
            {
                Matrix M = A;
                for (int i = 1; i < k; i++)
                {
                    M *= A;
                }

                return M;
            }
        }


        #endregion

        #region Virtual鰏

        /// <summary>
        /// Access the component in row i, column j of a non-empty matrix.
        /// </summary>
        /// <param name="i">One-based row index.</param>
        /// <param name="j">One-based column index.</param>
        /// <returns></returns>
        public virtual Complex this[int i, int j]
        {
            set
            {
                if (i <= 0 || j <= 0)
                    throw new ArgumentOutOfRangeException("Indices must be real positive.");

                if (i > rowCount)
                {
                    // dynamically add i-Rows new rows...
                    for (int k = 0; k < i - rowCount; k++)
                    {
                        this.Values.Add(new ArrayList(columnCount));

                        // ...with Cols columns
                        for (int t = 0; t < columnCount; t++)
                        {
                            ((ArrayList)Values[rowCount + k]).Add(Complex.Zero);
                        }
                    }

                    rowCount = i; // ha!
                }


                if (j > columnCount)
                {
                    // dynamically add j-Cols columns to each row
                    for (int k = 0; k < rowCount; k++)
                    {
                        for (int t = 0; t < j - columnCount; t++)
                        {
                            ((ArrayList)Values[k]).Add(Complex.Zero);
                        }
                    }

                    columnCount = j;
                }

                ((ArrayList)Values[i - 1])[j - 1] = value;
                //this.Values[i - 1, j - 1] = value; 
            }
            get
            {
                if (i > 0 && i <= rowCount && j > 0 && j <= columnCount)
                {
                    //Complex buf;
                    //
                    //try
                    //{
                    //    buf = (Complex)(((ArrayList)Values[i - 1])[j - 1]);
                    //}
                    //catch
                    //{
                    //    buf = new Complex((double)((int)(((ArrayList)Values[i - 1])[j - 1])));
                    //}
                    //
                    // return buf;                    

                    return (Complex)(((ArrayList)Values[i - 1])[j - 1]);
                }
                else
                    throw new ArgumentOutOfRangeException("Indices must not exceed size of matrix.");
            }
        }

        /// <summary>
        /// Access to the i-th component of an n by one matrix (column vector)
        /// or one by n matrix (row vector).
        /// </summary>
        /// <param name="i">One-based index.</param>
        /// <returns></returns>
        public virtual Complex this[int i]
        {
            set
            {
                if (rowCount == 1)
                {
                    // row vector

                    // dynamically extend vector if necessary
                    if (i > columnCount)
                    {
                        // dynamically add j-Cols columns to each row
                        for (int t = 0; t < i - columnCount; t++)
                            ((ArrayList)Values[0]).Add(Complex.Zero);

                        columnCount = i;
                    }

                    ((ArrayList)Values[0])[i - 1] = value;
                }
                else if (columnCount == 1)
                {
                    // column vector

                    if (i > rowCount)
                    {
                        // dynamically add i-Rows new rows...
                        for (int k = 0; k < i - rowCount; k++)
                        {
                            this.Values.Add(new ArrayList(columnCount));

                            // ...with one column each
                            ((ArrayList)Values[rowCount + k]).Add(Complex.Zero);
                        }

                        rowCount = i; // ha!
                    }

                    ((ArrayList)Values[i - 1])[0] = value;
                }
                else
                    throw new InvalidOperationException("Cannot access multidimensional matrix via single index.");
            }
            get
            {                                                    
                //Complex buf;

                //if (this.rowCount == 1)
                //    try
                //    {
                //        buf = (Complex)(((ArrayList)Values[0])[i - 1]);
                //    }
                //    catch
                //    {
                //        buf = new Complex((double)((int)(((ArrayList)Values[0])[i - 1])));
                //    }
                //else
                //    try
                //    {
                //        buf = (Complex)(((ArrayList)Values[i - 1])[0]);
                //    }
                //    catch
                //    {
                //        buf = new Complex((double)((int)(((ArrayList)Values[i - 1])[0])));
                //    }

                //return buf;


                if (this.RowCount == 1) // row vector
                    return (Complex)(((ArrayList)Values[0])[i - 1]);
                else if (this.ColumnCount == 1) // coumn vector
                    return (Complex)(((ArrayList)Values[i - 1])[0]);
                else // neither
                    throw new InvalidOperationException("General matrix acces requires double indexing.");
            }
        }
     
       
        #endregion

        #region Structs & Enums

        public enum DefinitenessType
        { 
            PositiveDefinite,
            PositiveSemidefinite,
            NegativeDefinite,
            NegativeSemidefinite,
            Indefinite
        }

        #endregion
    }
}

//namespace TF300.App.GUI.DatabaseUI.Utility
//{
//    /**
//    * 操作矩阵的类 Matrix

//    * @author 周长发
//    * @version 1.0
//    */
//    public class Matrix 
//    {
//        private int	numColumns = 0;			// 矩阵列数
//        private int	numRows = 0;			// 矩阵行数
//        private double eps = 0.0;			// 缺省精度
//        private double[] elements = null;	// 矩阵数据缓冲区

//        /**
//         * 属性: 矩阵列数
//         */
//        public int Columns
//        {
//            get
//            {
//                return numColumns;
//            }
//        }

//        /**
//         * 属性: 矩阵行数
//         */
//        public int Rows
//        {
//            get
//            {
//                return numRows;
//            }
//        }

//        /**
//         * 索引器: 访问矩阵元素
//         * @param row - 元素的行
//         * @param col - 元素的列
//         */
//        public double this[int row, int col]
//        {
//            get
//            {
//                return elements[col + row * numColumns]; 
//            }
//            set
//            {
//                elements[col + row * numColumns] = value;
//            }
//        }

//        /**
//         * 属性: Eps
//         */
//        public double Eps
//        {
//            get
//            {
//                return eps;
//            }
//            set
//            {
//                eps = value;
//            }
//        }

//        /**
//         * 基本构造函数
//         */
//        public Matrix()
//        {
//            numColumns = 1;
//            numRows = 1;
//            Init(numRows, numColumns);
//        }

//        /**
//         * 指定行列构造函数
//         * 
//         * @param nRows - 指定的矩阵行数
//         * @param nCols - 指定的矩阵列数
//         */
//        public Matrix(int nRows, int nCols)
//        {
//            numRows = nRows;
//            numColumns = nCols;
//            Init(numRows, numColumns);
//        }

//        /**
//         * 指定值构造函数
//         * 
//         * @param value - 二维数组，存储矩阵各元素的值
//         */
//        public Matrix(double[,] value)
//        {
//            numRows = value.GetLength(0);
//            numColumns = value.GetLength(1);
//            double[] data = new double[numRows * numColumns];
//            int k = 0;
//            for (int i=0; i<numRows; ++i)
//            {
//                for (int j=0; j<numColumns; ++j)
//                {
//                    data[k++] = value[i, j]; 
//                }
//            }
//            Init(numRows, numColumns);
//            SetData(data);
//        }

//        /**
//         * 指定值构造函数
//         * 
//         * @param nRows - 指定的矩阵行数
//         * @param nCols - 指定的矩阵列数
//         * @param value - 一维数组，长度为nRows*nCols，存储矩阵各元素的值
//         */
//        public Matrix(int nRows, int nCols, double[] value)
//        {
//            numRows = nRows;
//            numColumns = nCols;
//            Init(numRows, numColumns);
//            SetData(value);
//        }

//        /**
//         * 方阵构造函数
//         * 
//         * @param nSize - 方阵行列数
//         */
//        public Matrix(int nSize)
//        {
//            numRows = nSize;
//            numColumns = nSize;
//            Init(nSize, nSize);
//        }

//        /**
//         * 方阵构造函数
//         * 
//         * @param nSize - 方阵行列数
//         * @param value - 一维数组，长度为nRows*nRows，存储方阵各元素的值
//         */
//        public Matrix(int nSize, double[] value)
//        {
//            numRows = nSize;
//            numColumns = nSize;
//            Init(nSize, nSize);
//            SetData(value);
//        }

//        /**
//         * 拷贝构造函数
//         * 
//         * @param other - 源矩阵
//         */
//        public Matrix( Matrix other)
//        {
//            numColumns = other.GetNumColumns();
//            numRows = other.GetNumRows();
//            Init(numRows, numColumns);
//            SetData(other.elements);
//        }

//        /**
//         * 初始化函数
//         * 
//         * @param nRows - 指定的矩阵行数
//         * @param nCols - 指定的矩阵列数
//         * @return bool, 成功返回true, 否则返回false
//         */
//        public bool Init(int nRows, int nCols)
//        {
//            numRows = nRows;
//            numColumns = nCols;
//            int nSize = nCols*nRows;
//            if (nSize < 0)
//                return false;

//            // 分配内存
//            elements = new double[nSize];
		
//            return true;
//        }

//        /**
//         * 设置矩阵运算的精度
//         * 
//         * @param newEps - 新的精度值
//         */
//        public void SetEps(double newEps)
//        {
//            eps = newEps;
//        }
	
//        /**
//         * 取矩阵的精度值
//         * 
//         * @return double型，矩阵的精度值
//         */
//        public double GetEps()
//        {
//            return eps;
//        }

//        /**
//         * 重载 + 运算符
//         * 
//         * @return Matrix对象
//         */
//        public static Matrix operator +(Matrix m1, Matrix m2)
//        {
//            return m1.Add(m2);
//        }

//        /**
//         * 重载 - 运算符
//         * 
//         * @return Matrix对象
//         */
//        public static Matrix operator -(Matrix m1, Matrix m2)
//        {
//            return m1.Subtract(m2);
//        }

//        /**
//         * 重载 * 运算符
//         * 
//         * @return Matrix对象
//         */
//        public static Matrix operator *(Matrix m1, Matrix m2)
//        {
//            return m1.Multiply(m2);
//        }

//        /**
//         * 重载 double[] 运算符
//         * 
//         * @return double[]对象
//         */
//        public static implicit operator double[](Matrix m)
//        {
//            return m.elements;
//        }

//        /**
//         * 将方阵初始化为单位矩阵
//         * 
//         * @param nSize - 方阵行列数
//         * @return bool 型，初始化是否成功
//         */
//        public bool MakeUnitMatrix(int nSize)
//        {
//            if (! Init(nSize, nSize))
//                return false;

//            for (int i=0; i<nSize; ++i)
//                for (int j=0; j<nSize; ++j)
//                    if (i == j)
//                        SetElement(i, j, 1);

//            return true;
//        }

//        /**
//         * 将矩阵各元素的值转化为字符串, 元素之间的分隔符为",", 行与行之间有回车换行符
//         * @return string 型，转换得到的字符串
//         */
//        public override string ToString() 
//        {
//            return ToString(",", true);
//        }
	
//        /**
//         * 将矩阵各元素的值转化为字符串
//         * 
//         * @param sDelim - 元素之间的分隔符
//         * @param bLineBreak - 行与行之间是否有回车换行符
//         * @return string 型，转换得到的字符串
//         */
//        public string ToString(string sDelim, bool bLineBreak) 
//        {
//            string s="";

//            for (int i=0; i<numRows; ++i)
//            {
//                for (int j=0; j<numColumns; ++j)
//                {
//                    string ss = GetElement(i, j).ToString("F");
//                    s += ss;

//                    if (bLineBreak)
//                    {
//                        if (j != numColumns-1)
//                            s += sDelim;
//                    }
//                    else
//                    {
//                        if (i != numRows-1 || j != numColumns-1)
//                            s += sDelim;
//                    }
//                }
//                if (bLineBreak)
//                    if (i != numRows-1)
//                        s += "\r\n";
//            }

//            return s;
//        }

//        /**
//         * 将矩阵指定行中各元素的值转化为字符串
//         * 
//         * @param nRow - 指定的矩阵行，nRow = 0表示第一行
//         * @param sDelim - 元素之间的分隔符
//         * @return string 型，转换得到的字符串
//         */
//        public string ToStringRow(int nRow,  string sDelim) 
//        {
//            string s = "";

//            if (nRow >= numRows)
//                return s;

//            for (int j=0; j<numColumns; ++j)
//            {
//                string ss = GetElement(nRow, j).ToString("F");
//                s += ss;
//                if (j != numColumns-1)
//                    s += sDelim;
//            }

//            return s;
//        }

//        /**
//         * 将矩阵指定列中各元素的值转化为字符串
//         * 
//         * @param nCol - 指定的矩阵行，nCol = 0表示第一列
//         * @param sDelim - 元素之间的分隔符
//         * @return string 型，转换得到的字符串
//         */
//        public string ToStringCol(int nCol,  string sDelim /*= " "*/) 
//        {
//            string s = "";

//            if (nCol >= numColumns)
//                return s;

//            for (int i=0; i<numRows; ++i)
//            {
//                string ss = GetElement(i, nCol).ToString("F");
//                s += ss;
//                if (i != numRows-1)
//                    s += sDelim;
//            }

//            return s;
//        }

//        /**
//         * 设置矩阵各元素的值
//         * 
//         * @param value - 一维数组，长度为numColumns*numRows，存储
//         *	              矩阵各元素的值
//         */
//        public void SetData(double[] value)
//        {
//            elements = (double[])value.Clone();
//        }

//        /**
//         * 设置指定元素的值
//         * 
//         * @param nRow - 元素的行
//         * @param nCol - 元素的列
//         * @param value - 指定元素的值
//         * @return bool 型，说明设置是否成功
//         */
//        public bool SetElement(int nRow, int nCol, double value)
//        {
//            if (nCol < 0 || nCol >= numColumns || nRow < 0 || nRow >= numRows)
//                return false;						// array bounds error
		
//            elements[nCol + nRow * numColumns] = value;

//            return true;
//        }

//        /**
//         * 获取指定元素的值
//         * 
//         * @param nRow - 元素的行
//         * @param nCol - 元素的列
//         * @return double 型，指定元素的值
//         */
//        public double GetElement(int nRow, int nCol) 
//        {
//            return elements[nCol + nRow * numColumns] ;
//        }

//        /**
//         * 获取矩阵的列数
//         * 
//         * @return int 型，矩阵的列数
//         */
//        public int	GetNumColumns() 
//        {
//            return numColumns;
//        }

//        /**
//         * 获取矩阵的行数
//         * @return int 型，矩阵的行数
//         */
//        public int	GetNumRows() 
//        {
//            return numRows;
//        }

//        /**
//         * 获取矩阵的数据
//         * 
//         * @return double型数组，指向矩阵各元素的数据缓冲区
//         */
//        public double[] GetData() 
//        {
//            return elements;
//        }

//        /**
//         * 获取指定行的向量
//         * 
//         * @param nRow - 向量所在的行
//         * @param pVector - 指向向量中各元素的缓冲区
//         * @return int 型，向量中元素的个数，即矩阵的列数
//         */
//        public int GetRowVector(int nRow, double[] pVector) 
//        {
//            for (int j=0; j<numColumns; ++j)
//                pVector[j] = GetElement(nRow, j);

//            return numColumns;
//        }

//        /**
//         * 获取指定列的向量
//         * 
//         * @param nCol - 向量所在的列
//         * @param pVector - 指向向量中各元素的缓冲区
//         * @return int 型，向量中元素的个数，即矩阵的行数
//         */
//        public int GetColVector(int nCol, double[] pVector) 
//        {
//            for (int i=0; i<numRows; ++i)
//                pVector[i] = GetElement(i, nCol);

//            return numRows;
//        }

//        /**
//         * 给矩阵赋值
//         * 
//         * @param other - 用于给矩阵赋值的源矩阵
//         * @return Matrix型，阵与other相等
//         */
//        public Matrix SetValue(Matrix other)
//        {
//            if (other != this)
//            {
//                Init(other.GetNumRows(), other.GetNumColumns());
//                SetData(other.elements);
//            }

//            // finally return a reference to ourselves
//            return this ;
//        }

//        /**
//         * 判断矩阵否相等
//         * 
//         * @param other - 用于比较的矩阵
//         * @return bool 型，两个矩阵相等则为true，否则为false
//         */
//        public override bool Equals(object other) 
//        {
//            Matrix matrix = other as Matrix;
//            if (matrix == null)
//                return false;

//            // 首先检查行列数是否相等
//            if (numColumns != matrix.GetNumColumns() || numRows != matrix.GetNumRows())
//                return false;

//            for (int i=0; i<numRows; ++i)
//            {
//                for (int j=0; j<numColumns; ++j)
//                {
//                    if (Math.Abs(GetElement(i, j) - matrix.GetElement(i, j)) > eps)
//                        return false;
//                }
//            }

//            return true;
//        }

//        /**
//         * 因为重写了Equals，因此必须重写GetHashCode
//         * 
//         * @return int型，返回复数对象散列码
//         */
//        public override int GetHashCode()
//        {
//            double sum = 0;
//            for (int i=0; i<numRows; ++i)
//            {
//                for (int j=0; j<numColumns; ++j)
//                {
//                    sum += Math.Abs(GetElement(i, j));
//                }
//            }
//            return (int)Math.Sqrt(sum);
//        }

//        /**
//         * 实现矩阵的加法
//         * 
//         * @param other - 与指定矩阵相加的矩阵
//         * @return Matrix型，指定矩阵与other相加之和
//         * @如果矩阵的行/列数不匹配，则会抛出异常
//         */
//        public Matrix Add(Matrix other) 
//        {
//            // 首先检查行列数是否相等
//            if (numColumns != other.GetNumColumns() ||
//                numRows != other.GetNumRows())
//                throw new Exception("矩阵的行/列数不匹配。");

//            // 构造结果矩阵
//            Matrix	result = new Matrix(this) ;		// 拷贝构造
		
//            // 矩阵加法
//            for (int i = 0 ; i < numRows ; ++i)
//            {
//                for (int j = 0 ; j <  numColumns; ++j)
//                    result.SetElement(i, j, result.GetElement(i, j) + other.GetElement(i, j)) ;
//            }

//            return result ;
//        }

//        /**
//         * 实现矩阵的减法
//         * 
//         * @param other - 与指定矩阵相减的矩阵
//         * @return Matrix型，指定矩阵与other相减之差
//         * @如果矩阵的行/列数不匹配，则会抛出异常
//         */
//        public Matrix Subtract(Matrix other) 
//        {
//            if (numColumns != other.GetNumColumns() ||
//                numRows != other.GetNumRows())
//                throw new Exception("矩阵的行/列数不匹配。");

//            // 构造结果矩阵
//            Matrix	result = new Matrix(this) ;		// 拷贝构造

//            // 进行减法操作
//            for (int i = 0 ; i < numRows ; ++i)
//            {
//                for (int j = 0 ; j <  numColumns; ++j)
//                    result.SetElement(i, j, result.GetElement(i, j) - other.GetElement(i, j)) ;
//            }

//            return result ;
//        }

//        /**
//         * 实现矩阵的数乘
//         * 
//         * @param value - 与指定矩阵相乘的实数
//         * @return Matrix型，指定矩阵与value相乘之积
//         */
//        public Matrix Multiply(double value) 
//        {
//            // 构造目标矩阵
//            Matrix	result = new Matrix(this) ;		// copy ourselves
		
//            // 进行数乘
//            for (int i = 0 ; i < numRows ; ++i)
//            {
//                for (int j = 0 ; j <  numColumns; ++j)
//                    result.SetElement(i, j, result.GetElement(i, j) * value) ;
//            }

//            return result ;
//        }

//        /**
//         * 实现矩阵的乘法
//         * 
//         * @param other - 与指定矩阵相乘的矩阵
//         * @return Matrix型，指定矩阵与other相乘之积
//         * @如果矩阵的行/列数不匹配，则会抛出异常
//         */
//        public Matrix Multiply(Matrix other) 
//        {
//            // 首先检查行列数是否符合要求
//            if (numColumns != other.GetNumRows())
//                throw new Exception("矩阵的行/列数不匹配。");

//            // ruct the object we are going to return
//            Matrix	result = new Matrix(numRows, other.GetNumColumns());

//            // 矩阵乘法，即
//            //
//            // [A][B][C]   [G][H]     [A*G + B*I + C*K][A*H + B*J + C*L]
//            // [D][E][F] * [I][J] =   [D*G + E*I + F*K][D*H + E*J + F*L]
//            //             [K][L]
//            //
//            double	value ;
//            for (int i = 0 ; i < result.GetNumRows() ; ++i)
//            {
//                for (int j = 0 ; j < other.GetNumColumns() ; ++j)
//                {
//                    value = 0.0 ;
//                    for (int k = 0 ; k < numColumns ; ++k)
//                    {
//                        value += GetElement(i, k) * other.GetElement(k, j) ;
//                    }

//                    result.SetElement(i, j, value) ;
//                }
//            }

//            return result ;
//        }

//        /**
//         * 复矩阵的乘法
//         * 
//         * @param AR - 左边复矩阵的实部矩阵
//         * @param AI - 左边复矩阵的虚部矩阵
//         * @param BR - 右边复矩阵的实部矩阵
//         * @param BI - 右边复矩阵的虚部矩阵
//         * @param CR - 乘积复矩阵的实部矩阵
//         * @param CI - 乘积复矩阵的虚部矩阵
//         * @return bool型，复矩阵乘法是否成功
//         */
//        public bool Multiply(Matrix AR,  Matrix AI,  Matrix BR,  Matrix BI, Matrix CR, Matrix CI) 
//        {
//            // 首先检查行列数是否符合要求
//            if (AR.GetNumColumns() != AI.GetNumColumns() ||
//                AR.GetNumRows() != AI.GetNumRows() ||
//                BR.GetNumColumns() != BI.GetNumColumns() ||
//                BR.GetNumRows() != BI.GetNumRows() ||
//                AR.GetNumColumns() != BR.GetNumRows())
//                return false;

//            // 构造乘积矩阵实部矩阵和虚部矩阵
//            Matrix mtxCR = new Matrix(AR.GetNumRows(), BR.GetNumColumns());
//            Matrix mtxCI = new Matrix(AR.GetNumRows(), BR.GetNumColumns());
//            // 复矩阵相乘
//            for (int i=0; i<AR.GetNumRows(); ++i)
//            {
//                for (int j=0; j<BR.GetNumColumns(); ++j)
//                {
//                    double vr = 0;
//                    double vi = 0;
//                    for (int k =0; k<AR.GetNumColumns(); ++k)
//                    {
//                        double p = AR.GetElement(i, k) * BR.GetElement(k, j);
//                        double q = AI.GetElement(i, k) * BI.GetElement(k, j);
//                        double s = (AR.GetElement(i, k) + AI.GetElement(i, k)) * (BR.GetElement(k, j) + BI.GetElement(k, j));
//                        vr += p - q;
//                        vi += s - p - q;
//                    }
//                    mtxCR.SetElement(i, j, vr);
//                    mtxCI.SetElement(i, j, vi);
//                }
//            }

//            CR = mtxCR;
//            CI = mtxCI;

//            return true;
//        }

//        /**
//         * 矩阵的转置
//         * 
//         * @return Matrix型，指定矩阵转置矩阵
//         */
//        public Matrix Transpose() 
//        {
//            // 构造目标矩阵
//            Matrix	Trans = new Matrix(numColumns, numRows);

//            // 转置各元素
//            for (int i = 0 ; i < numRows ; ++i)
//            {
//                for (int j = 0 ; j < numColumns ; ++j)
//                    Trans.SetElement(j, i, GetElement(i, j)) ;
//            }

//            return Trans;
//        }

//        /**
//         * 实矩阵求逆的全选主元高斯－约当法
//         * 
//         * @return bool型，求逆是否成功
//         */
//        public bool InvertGaussJordan()
//        {
//            int i,j,k,l,u,v;
//            double d = 0, p = 0;

//            // 分配内存
//            int[] pnRow = new int[numColumns];
//            int[] pnCol = new int[numColumns];

//            // 消元
//            for (k=0; k<=numColumns-1; k++)
//            { 
//                d=0.0;
//                for (i=k; i<=numColumns-1; i++)
//                {
//                    for (j=k; j<=numColumns-1; j++)
//                    { 
//                        l=i*numColumns+j; p=Math.Abs(elements[l]);
//                        if (p>d) 
//                        { 
//                            d=p; 
//                            pnRow[k]=i; 
//                            pnCol[k]=j;
//                        }
//                    }
//                }
	        
//                // 失败
//                if (d == 0.0)
//                {
//                    return false;
//                }

//                if (pnRow[k] != k)
//                {
//                    for (j=0; j<=numColumns-1; j++)
//                    { 
//                        u=k*numColumns+j; 
//                        v=pnRow[k]*numColumns+j;
//                        p=elements[u]; 
//                        elements[u]=elements[v]; 
//                        elements[v]=p;
//                    }
//                }
	        
//                if (pnCol[k] != k)
//                {
//                    for (i=0; i<=numColumns-1; i++)
//                    { 
//                        u=i*numColumns+k; 
//                        v=i*numColumns+pnCol[k];
//                        p=elements[u]; 
//                        elements[u]=elements[v]; 
//                        elements[v]=p;
//                    }
//                }

//                l=k*numColumns+k;
//                elements[l]=1.0/elements[l];
//                for (j=0; j<=numColumns-1; j++)
//                {
//                    if (j != k)
//                    { 
//                        u=k*numColumns+j; 
//                        elements[u]=elements[u]*elements[l];
//                    }
//                }

//                for (i=0; i<=numColumns-1; i++)
//                {
//                    if (i!=k)
//                    {
//                        for (j=0; j<=numColumns-1; j++)
//                        {
//                            if (j!=k)
//                            { 
//                                u=i*numColumns+j;
//                                elements[u]=elements[u]-elements[i*numColumns+k]*elements[k*numColumns+j];
//                            }
//                        }
//                    }
//                }

//                for (i=0; i<=numColumns-1; i++)
//                {
//                    if (i!=k)
//                    { 
//                        u=i*numColumns+k; 
//                        elements[u]=-elements[u]*elements[l];
//                    }
//                }
//            }

//            // 调整恢复行列次序
//            for (k=numColumns-1; k>=0; k--)
//            { 
//                if (pnCol[k]!=k)
//                {
//                    for (j=0; j<=numColumns-1; j++)
//                    { 
//                        u=k*numColumns+j; 
//                        v=pnCol[k]*numColumns+j;
//                        p=elements[u]; 
//                        elements[u]=elements[v]; 
//                        elements[v]=p;
//                    }
//                }

//                if (pnRow[k]!=k)
//                {
//                    for (i=0; i<=numColumns-1; i++)
//                    { 
//                        u=i*numColumns+k; 
//                        v=i*numColumns+pnRow[k];
//                        p=elements[u]; 
//                        elements[u]=elements[v]; 
//                        elements[v]=p;
//                    }
//                }
//            }

//            // 成功返回
//            return true;
//        }

//        /**
//         * 复矩阵求逆的全选主元高斯－约当法
//         * 
//         * @param mtxImag - 复矩阵的虚部矩阵，当前矩阵为复矩阵的实部
//         * @return bool型，求逆是否成功
//         */
//        public bool InvertGaussJordan(Matrix mtxImag)
//        {
//            int i,j,k,l,u,v,w;
//            double p,q,s,t,d,b;

//            // 分配内存
//            int[] pnRow = new int[numColumns];
//            int[] pnCol = new int[numColumns];

//            // 消元
//            for (k=0; k<=numColumns-1; k++)
//            { 
//                d=0.0;
//                for (i=k; i<=numColumns-1; i++)
//                {
//                    for (j=k; j<=numColumns-1; j++)
//                    { 
//                        u=i*numColumns+j;
//                        p=elements[u]*elements[u]+mtxImag.elements[u]*mtxImag.elements[u];
//                        if (p>d) 
//                        { 
//                            d=p; 
//                            pnRow[k]=i; 
//                            pnCol[k]=j;
//                        }
//                    }
//                }

//                // 失败
//                if (d == 0.0)
//                { 
//                    return false;
//                }

//                if (pnRow[k]!=k)
//                {
//                    for (j=0; j<=numColumns-1; j++)
//                    { 
//                        u=k*numColumns+j; 
//                        v=pnRow[k]*numColumns+j;
//                        t=elements[u]; 
//                        elements[u]=elements[v]; 
//                        elements[v]=t;
//                        t=mtxImag.elements[u]; 
//                        mtxImag.elements[u]=mtxImag.elements[v]; 
//                        mtxImag.elements[v]=t;
//                    }
//                }

//                if (pnCol[k]!=k)
//                {
//                    for (i=0; i<=numColumns-1; i++)
//                    { 
//                        u=i*numColumns+k; 
//                        v=i*numColumns+pnCol[k];
//                        t=elements[u]; 
//                        elements[u]=elements[v]; 
//                        elements[v]=t;
//                        t=mtxImag.elements[u]; 
//                        mtxImag.elements[u]=mtxImag.elements[v]; 
//                        mtxImag.elements[v]=t;
//                    }
//                }

//                l=k*numColumns+k;
//                elements[l]=elements[l]/d; mtxImag.elements[l]=-mtxImag.elements[l]/d;
//                for (j=0; j<=numColumns-1; j++)
//                {
//                    if (j!=k)
//                    { 
//                        u=k*numColumns+j;
//                        p=elements[u]*elements[l]; 
//                        q=mtxImag.elements[u]*mtxImag.elements[l];
//                        s=(elements[u]+mtxImag.elements[u])*(elements[l]+mtxImag.elements[l]);
//                        elements[u]=p-q; 
//                        mtxImag.elements[u]=s-p-q;
//                    }
//                }

//                for (i=0; i<=numColumns-1; i++)
//                {
//                    if (i!=k)
//                    { 
//                        v=i*numColumns+k;
//                        for (j=0; j<=numColumns-1; j++)
//                        {
//                            if (j!=k)
//                            { 
//                                u=k*numColumns+j;  
//                                w=i*numColumns+j;
//                                p=elements[u]*elements[v]; 
//                                q=mtxImag.elements[u]*mtxImag.elements[v];
//                                s=(elements[u]+mtxImag.elements[u])*(elements[v]+mtxImag.elements[v]);
//                                t=p-q; 
//                                b=s-p-q;
//                                elements[w]=elements[w]-t;
//                                mtxImag.elements[w]=mtxImag.elements[w]-b;
//                            }
//                        }
//                    }
//                }

//                for (i=0; i<=numColumns-1; i++)
//                {
//                    if (i!=k)
//                    { 
//                        u=i*numColumns+k;
//                        p=elements[u]*elements[l]; 
//                        q=mtxImag.elements[u]*mtxImag.elements[l];
//                        s=(elements[u]+mtxImag.elements[u])*(elements[l]+mtxImag.elements[l]);
//                        elements[u]=q-p; 
//                        mtxImag.elements[u]=p+q-s;
//                    }
//                }
//            }

//            // 调整恢复行列次序
//            for (k=numColumns-1; k>=0; k--)
//            { 
//                if (pnCol[k]!=k)
//                {
//                    for (j=0; j<=numColumns-1; j++)
//                    { 
//                        u=k*numColumns+j; 
//                        v=pnCol[k]*numColumns+j;
//                        t=elements[u]; 
//                        elements[u]=elements[v]; 
//                        elements[v]=t;
//                        t=mtxImag.elements[u]; 
//                        mtxImag.elements[u]=mtxImag.elements[v]; 
//                        mtxImag.elements[v]=t;
//                    }
//                }

//                if (pnRow[k]!=k)
//                {
//                    for (i=0; i<=numColumns-1; i++)
//                    { 
//                        u=i*numColumns+k; 
//                        v=i*numColumns+pnRow[k];
//                        t=elements[u]; 
//                        elements[u]=elements[v]; 
//                        elements[v]=t;
//                        t=mtxImag.elements[u]; 
//                        mtxImag.elements[u]=mtxImag.elements[v]; 
//                        mtxImag.elements[v]=t;
//                    }
//                }
//            }

//            // 成功返回
//            return true;
//        }

//        /**
//         * 对称正定矩阵的求逆
//         * 
//         * @return bool型，求逆是否成功
//         */
//        public bool InvertSsgj()
//        { 
//            int i, j ,k, m;
//            double w, g;

//            // 临时内存
//            double[] pTmp = new double[numColumns];

//            // 逐列处理
//            for (k=0; k<=numColumns-1; k++)
//            { 
//                w=elements[0];
//                if (w == 0.0)
//                { 
//                    return false;
//                }

//                m=numColumns-k-1;
//                for (i=1; i<=numColumns-1; i++)
//                { 
//                    g=elements[i*numColumns]; 
//                    pTmp[i]=g/w;
//                    if (i<=m) 
//                        pTmp[i]=-pTmp[i];
//                    for (j=1; j<=i; j++)
//                        elements[(i-1)*numColumns+j-1]=elements[i*numColumns+j]+g*pTmp[j];
//                }

//                elements[numColumns*numColumns-1]=1.0/w;
//                for (i=1; i<=numColumns-1; i++)
//                    elements[(numColumns-1)*numColumns+i-1]=pTmp[i];
//            }

//            // 行列调整
//            for (i=0; i<=numColumns-2; i++)
//                for (j=i+1; j<=numColumns-1; j++)
//                    elements[i*numColumns+j]=elements[j*numColumns+i];

//            return true;
//        }

//        /**
//         * 托伯利兹矩阵求逆的埃兰特方法
//         * 
//         * @return bool型，求逆是否成功
//         */
//        public bool InvertTrench()
//        { 
//            int i,j,k;
//            double a,s;

//            // 上三角元素
//            double[] t = new double[numColumns];
//            // 下三角元素
//            double[] tt = new double[numColumns];

//            // 上、下三角元素赋值
//            for (i=0; i<numColumns; ++i)
//            {
//                t[i] = GetElement(0, i);
//                tt[i] = GetElement(i, 0);
//            }

//            // 临时缓冲区
//            double[] c = new double[numColumns];
//            double[] r = new double[numColumns];
//            double[] p = new double[numColumns];

//            // 非Toeplitz矩阵，返回
//            if (t[0] == 0.0)
//            { 
//                return false;
//            }

//            a=t[0]; 
//            c[0]=tt[1]/t[0]; 
//            r[0]=t[1]/t[0];

//            for (k=0; k<=numColumns-3; k++)
//            { 
//                s=0.0;
//                for (j=1; j<=k+1; j++)
//                    s=s+c[k+1-j]*tt[j];

//                s=(s-tt[k+2])/a;
//                for (i=0; i<=k; i++)
//                    p[i]=c[i]+s*r[k-i];

//                c[k+1]=-s;
//                s=0.0;
//                for (j=1; j<=k+1; j++)
//                    s=s+r[k+1-j]*t[j];
	        
//                s=(s-t[k+2])/a;
//                for (i=0; i<=k; i++)
//                { 
//                    r[i]=r[i]+s*c[k-i];
//                    c[k-i]=p[k-i];
//                }

//                r[k+1]=-s;
//                a=0.0;
//                for (j=1; j<=k+2; j++)
//                    a=a+t[j]*c[j-1];

//                a=t[0]-a;

//                // 求解失败
//                if (a == 0.0)
//                { 
//                    return false;
//                }
//            }

//            elements[0]=1.0/a;
//            for (i=0; i<=numColumns-2; i++)
//            { 
//                k=i+1; 
//                j=(i+1)*numColumns;
//                elements[k]=-r[i]/a; 
//                elements[j]=-c[i]/a;
//            }

//            for (i=0; i<=numColumns-2; i++)
//            {
//                for (j=0; j<=numColumns-2; j++)
//                { 
//                    k=(i+1)*numColumns+j+1;
//                    elements[k]=elements[i*numColumns+j]-c[i]*elements[j+1];
//                    elements[k]=elements[k]+c[numColumns-j-2]*elements[numColumns-i-1];
//                }
//            }

//            return true;
//        }

//        /**
//         * 求行列式值的全选主元高斯消去法
//         * 
//         * @return double型，行列式的值
//         */
//        public double ComputeDetGauss()
//        { 
//            int i,j,k,nis = 0,js = 0,l,u,v;
//            double f,det,q,d;
	    
//            // 初值
//            f=1.0; 
//            det=1.0;
	    
//            // 消元
//            for (k=0; k<=numColumns-2; k++)
//            { 
//                q=0.0;
//                for (i=k; i<=numColumns-1; i++)
//                {
//                    for (j=k; j<=numColumns-1; j++)
//                    { 
//                        l=i*numColumns+j; 
//                        d=Math.Abs(elements[l]);
//                        if (d>q) 
//                        { 
//                            q=d; 
//                            nis=i; 
//                            js=j;
//                        }
//                    }
//                }

//                if (q == 0.0)
//                { 
//                    det=0.0; 
//                    return(det);
//                }
	        
//                if (nis!=k)
//                { 
//                    f=-f;
//                    for (j=k; j<=numColumns-1; j++)
//                    { 
//                        u=k*numColumns+j; 
//                        v=nis*numColumns+j;
//                        d=elements[u]; 
//                        elements[u]=elements[v]; 
//                        elements[v]=d;
//                    }
//                }
	        
//                if (js!=k)
//                { 
//                    f=-f;
//                    for (i=k; i<=numColumns-1; i++)
//                    {
//                        u=i*numColumns+js; 
//                        v=i*numColumns+k;
//                        d=elements[u]; 
//                        elements[u]=elements[v]; 
//                        elements[v]=d;
//                    }
//                }

//                l=k*numColumns+k;
//                det=det*elements[l];
//                for (i=k+1; i<=numColumns-1; i++)
//                { 
//                    d=elements[i*numColumns+k]/elements[l];
//                    for (j=k+1; j<=numColumns-1; j++)
//                    { 
//                        u=i*numColumns+j;
//                        elements[u]=elements[u]-d*elements[k*numColumns+j];
//                    }
//                }
//            }
	    
//            // 求值
//            det=f*det*elements[numColumns*numColumns-1];

//            return(det);
//        }

//        /**
//         * 求矩阵秩的全选主元高斯消去法
//         * 
//         * @return int型，矩阵的秩
//         */
//        public int ComputeRankGauss()
//        { 
//            int i,j,k,nn,nis = 0,js = 0,l,ll,u,v;
//            double q,d;
	    
//            // 秩小于等于行列数
//            nn = numRows;
//            if (numRows >= numColumns) 
//                nn = numColumns;

//            k=0;

//            // 消元求解
//            for (l=0; l<=nn-1; l++)
//            { 
//                q=0.0;
//                for (i=l; i<=numRows-1; i++)
//                {
//                    for (j=l; j<=numColumns-1; j++)
//                    { 
//                        ll=i*numColumns+j; 
//                        d=Math.Abs(elements[ll]);
//                        if (d>q) 
//                        { 
//                            q=d; 
//                            nis=i; 
//                            js=j;
//                        }
//                    }
//                }

//                if (q == 0.0) 
//                    return(k);

//                k=k+1;
//                if (nis!=l)
//                { 
//                    for (j=l; j<=numColumns-1; j++)
//                    { 
//                        u=l*numColumns+j; 
//                        v=nis*numColumns+j;
//                        d=elements[u]; 
//                        elements[u]=elements[v]; 
//                        elements[v]=d;
//                    }
//                }
//                if (js!=l)
//                { 
//                    for (i=l; i<=numRows-1; i++)
//                    { 
//                        u=i*numColumns+js; 
//                        v=i*numColumns+l;
//                        d=elements[u]; 
//                        elements[u]=elements[v]; 
//                        elements[v]=d;
//                    }
//                }
	        
//                ll=l*numColumns+l;
//                for (i=l+1; i<=numColumns-1; i++)
//                { 
//                    d=elements[i*numColumns+l]/elements[ll];
//                    for (j=l+1; j<=numColumns-1; j++)
//                    { 
//                        u=i*numColumns+j;
//                        elements[u]=elements[u]-d*elements[l*numColumns+j];
//                    }
//                }
//            }
	    
//            return(k);
//        }

//        /**
//         * 对称正定矩阵的乔里斯基分解与行列式的求值
//         * 
//         * @param realDetValue - 返回行列式的值
//         * @return bool型，求解是否成功
//         */
//        public bool ComputeDetCholesky(ref double realDetValue)
//        { 
//            int i,j,k,u,l;
//            double d;
	    
//            // 不满足求解要求
//            if (elements[0] <= 0.0)
//                return false;

//            // 乔里斯基分解

//            elements[0]=Math.Sqrt(elements[0]);
//            d=elements[0];

//            for (i=1; i<=numColumns-1; i++)
//            { 
//                u=i*numColumns; 
//                elements[u]=elements[u]/elements[0];
//            }
	    
//            for (j=1; j<=numColumns-1; j++)
//            { 
//                l=j*numColumns+j;
//                for (k=0; k<=j-1; k++)
//                { 
//                    u=j*numColumns+k; 
//                    elements[l]=elements[l]-elements[u]*elements[u];
//                }
	        
//                if (elements[l] <= 0.0)
//                    return false;

//                elements[l]=Math.Sqrt(elements[l]);
//                d=d*elements[l];
	        
//                for (i=j+1; i<=numColumns-1; i++)
//                { 
//                    u=i*numColumns+j;
//                    for (k=0; k<=j-1; k++)
//                        elements[u]=elements[u]-elements[i*numColumns+k]*elements[j*numColumns+k];
	            
//                    elements[u]=elements[u]/elements[l];
//                }
//            }
	    
//            // 行列式求值
//            realDetValue=d*d;
		
//            // 下三角矩阵
//            for (i=0; i<=numColumns-2; i++)
//                for (j=i+1; j<=numColumns-1; j++)
//                    elements[i*numColumns+j]=0.0;

//            return true;
//        }

//        /**
//         * 矩阵的三角分解，分解成功后，原矩阵将成为Q矩阵
//         * 
//         * @param mtxL - 返回分解后的L矩阵
//         * @param mtxU - 返回分解后的U矩阵
//         * @return bool型，求解是否成功
//         */
//        public bool SplitLU(Matrix mtxL, Matrix mtxU)
//        { 
//            int i,j,k,w,v,ll;
	    
//            // 初始化结果矩阵
//            if (! mtxL.Init(numColumns, numColumns) ||
//                ! mtxU.Init(numColumns, numColumns))
//                return false;

//            for (k=0; k<=numColumns-2; k++)
//            { 
//                ll=k*numColumns+k;
//                if (elements[ll] == 0.0)
//                    return false;

//                for (i=k+1; i<=numColumns-1; i++)
//                { 
//                    w=i*numColumns+k; 
//                    elements[w]=elements[w]/elements[ll];
//                }

//                for (i=k+1; i<=numColumns-1; i++)
//                { 
//                    w=i*numColumns+k;
//                    for (j=k+1; j<=numColumns-1; j++)
//                    { 
//                        v=i*numColumns+j;
//                        elements[v]=elements[v]-elements[w]*elements[k*numColumns+j];
//                    }
//                }
//            }
	    
//            for (i=0; i<=numColumns-1; i++)
//            {
//                for (j=0; j<i; j++)
//                { 
//                    w=i*numColumns+j; 
//                    mtxL.elements[w]=elements[w]; 
//                    mtxU.elements[w]=0.0;
//                }

//                w=i*numColumns+i;
//                mtxL.elements[w]=1.0; 
//                mtxU.elements[w]=elements[w];
	        
//                for (j=i+1; j<=numColumns-1; j++)
//                { 
//                    w=i*numColumns+j; 
//                    mtxL.elements[w]=0.0; 
//                    mtxU.elements[w]=elements[w];
//                }
//            }

//            return true;
//        }

//        /**
//         * 一般实矩阵的QR分解，分解成功后，原矩阵将成为R矩阵
//         * 
//         * @param mtxQ - 返回分解后的Q矩阵
//         * @return bool型，求解是否成功
//         */
//        public bool SplitQR(Matrix mtxQ)
//        { 
//            int i,j,k,l,nn,p,jj;
//            double u,alpha,w,t;
	    
//            if (numRows < numColumns)
//                return false;

//            // 初始化Q矩阵
//            if (! mtxQ.Init(numRows, numRows))
//                return false;

//            // 对角线元素单位化
//            for (i=0; i<=numRows-1; i++)
//            {
//                for (j=0; j<=numRows-1; j++)
//                { 
//                    l=i*numRows+j; 
//                    mtxQ.elements[l]=0.0;
//                    if (i==j) 
//                        mtxQ.elements[l]=1.0;
//                }
//            }

//            // 开始分解

//            nn=numColumns;
//            if (numRows == numColumns) 
//                nn=numRows-1;

//            for (k=0; k<=nn-1; k++)
//            { 
//                u=0.0; 
//                l=k*numColumns+k;
//                for (i=k; i<=numRows-1; i++)
//                { 
//                    w=Math.Abs(elements[i*numColumns+k]);
//                    if (w>u) 
//                        u=w;
//                }
	        
//                alpha=0.0;
//                for (i=k; i<=numRows-1; i++)
//                { 
//                    t=elements[i*numColumns+k]/u; 
//                    alpha=alpha+t*t;
//                }

//                if (elements[l]>0.0) 
//                    u=-u;

//                alpha=u*Math.Sqrt(alpha);
//                if (alpha == 0.0)
//                    return false;

//                u=Math.Sqrt(2.0*alpha*(alpha-elements[l]));
//                if ((u+1.0)!=1.0)
//                { 
//                    elements[l]=(elements[l]-alpha)/u;
//                    for (i=k+1; i<=numRows-1; i++)
//                    { 
//                        p=i*numColumns+k; 
//                        elements[p]=elements[p]/u;
//                    }
	            
//                    for (j=0; j<=numRows-1; j++)
//                    { 
//                        t=0.0;
//                        for (jj=k; jj<=numRows-1; jj++)
//                            t=t+elements[jj*numColumns+k]*mtxQ.elements[jj*numRows+j];

//                        for (i=k; i<=numRows-1; i++)
//                        { 
//                            p=i*numRows+j; 
//                            mtxQ.elements[p]=mtxQ.elements[p]-2.0*t*elements[i*numColumns+k];
//                        }
//                    }
	            
//                    for (j=k+1; j<=numColumns-1; j++)
//                    { 
//                        t=0.0;
	                
//                        for (jj=k; jj<=numRows-1; jj++)
//                            t=t+elements[jj*numColumns+k]*elements[jj*numColumns+j];
	                
//                        for (i=k; i<=numRows-1; i++)
//                        { 
//                            p=i*numColumns+j; 
//                            elements[p]=elements[p]-2.0*t*elements[i*numColumns+k];
//                        }
//                    }
	            
//                    elements[l]=alpha;
//                    for (i=k+1; i<=numRows-1; i++)
//                        elements[i*numColumns+k]=0.0;
//                }
//            }
	    
//            // 调整元素
//            for (i=0; i<=numRows-2; i++)
//            {
//                for (j=i+1; j<=numRows-1;j++)
//                { 
//                    p=i*numRows+j; 
//                    l=j*numRows+i;
//                    t=mtxQ.elements[p]; 
//                    mtxQ.elements[p]=mtxQ.elements[l]; 
//                    mtxQ.elements[l]=t;
//                }
//            }

//            return true;
//        }

//        /**
//         * 一般实矩阵的奇异值分解，分解成功后，原矩阵对角线元素就是矩阵的奇异值
//         * 
//         * @param mtxU - 返回分解后的U矩阵
//         * @param mtxV - 返回分解后的V矩阵
//         * @param eps - 计算精度
//         * @return bool型，求解是否成功
//         */
//        public bool SplitUV(Matrix mtxU, Matrix mtxV, double eps)
//        { 
//            int i,j,k,l,it,ll,kk,ix,iy,mm,nn,iz,m1,ks;
//            double d,dd,t,sm,sm1,em1,sk,ek,b,c,shh;
//            double[] fg = new double[2];
//            double[] cs = new double[2];

//            int m = numRows;
//            int n = numColumns;

//            // 初始化U, V矩阵
//            if (! mtxU.Init(m, m) || ! mtxV.Init(n, n))
//                return false;

//            // 临时缓冲区
//            int ka = Math.Max(m, n) + 1;
//            double[] s = new double[ka];
//            double[] e = new double[ka];
//            double[] w = new double[ka];

//            // 指定迭代次数为60
//            it=60; 
//            k=n;

//            if (m-1<n) 
//                k=m-1;

//            l=m;
//            if (n-2<m) 
//                l=n-2;
//            if (l<0) 
//                l=0;

//            // 循环迭代计算
//            ll=k;
//            if (l>k) 
//                ll=l;
//            if (ll>=1)
//            { 
//                for (kk=1; kk<=ll; kk++)
//                { 
//                    if (kk<=k)
//                    { 
//                        d=0.0;
//                        for (i=kk; i<=m; i++)
//                        { 
//                            ix=(i-1)*n+kk-1; 
//                            d=d+elements[ix]*elements[ix];
//                        }

//                        s[kk-1]=Math.Sqrt(d);
//                        if (s[kk-1]!=0.0)
//                        { 
//                            ix=(kk-1)*n+kk-1;
//                            if (elements[ix]!=0.0)
//                            { 
//                                s[kk-1]=Math.Abs(s[kk-1]);
//                                if (elements[ix]<0.0) 
//                                    s[kk-1]=-s[kk-1];
//                            }
	                    
//                            for (i=kk; i<=m; i++)
//                            { 
//                                iy=(i-1)*n+kk-1;
//                                elements[iy]=elements[iy]/s[kk-1];
//                            }
	                    
//                            elements[ix]=1.0+elements[ix];
//                        }
	                
//                        s[kk-1]=-s[kk-1];
//                    }
	            
//                    if (n>=kk+1)
//                    { 
//                        for (j=kk+1; j<=n; j++)
//                        { 
//                            if ((kk<=k)&&(s[kk-1]!=0.0))
//                            { 
//                                d=0.0;
//                                for (i=kk; i<=m; i++)
//                                { 
//                                    ix=(i-1)*n+kk-1;
//                                    iy=(i-1)*n+j-1;
//                                    d=d+elements[ix]*elements[iy];
//                                }
	                        
//                                d=-d/elements[(kk-1)*n+kk-1];
//                                for (i=kk; i<=m; i++)
//                                { 
//                                    ix=(i-1)*n+j-1;
//                                    iy=(i-1)*n+kk-1;
//                                    elements[ix]=elements[ix]+d*elements[iy];
//                                }
//                            }
	                    
//                            e[j-1]=elements[(kk-1)*n+j-1];
//                        }
//                    }
	            
//                    if (kk<=k)
//                    { 
//                        for (i=kk; i<=m; i++)
//                        { 
//                            ix=(i-1)*m+kk-1; 
//                            iy=(i-1)*n+kk-1;
//                            mtxU.elements[ix]=elements[iy];
//                        }
//                    }
	            
//                    if (kk<=l)
//                    { 
//                        d=0.0;
//                        for (i=kk+1; i<=n; i++)
//                            d=d+e[i-1]*e[i-1];
	                
//                        e[kk-1]=Math.Sqrt(d);
//                        if (e[kk-1]!=0.0)
//                        { 
//                            if (e[kk]!=0.0)
//                            { 
//                                e[kk-1]=Math.Abs(e[kk-1]);
//                                if (e[kk]<0.0) 
//                                    e[kk-1]=-e[kk-1];
//                            }

//                            for (i=kk+1; i<=n; i++)
//                                e[i-1]=e[i-1]/e[kk-1];
	                    
//                            e[kk]=1.0+e[kk];
//                        }
	                
//                        e[kk-1]=-e[kk-1];
//                        if ((kk+1<=m)&& (e[kk-1]!=0.0))
//                        { 
//                            for (i=kk+1; i<=m; i++) 
//                                w[i-1]=0.0;
	                    
//                            for (j=kk+1; j<=n; j++)
//                                for (i=kk+1; i<=m; i++)
//                                    w[i-1]=w[i-1]+e[j-1]*elements[(i-1)*n+j-1];
	                    
//                            for (j=kk+1; j<=n; j++)
//                            {
//                                for (i=kk+1; i<=m; i++)
//                                { 
//                                    ix=(i-1)*n+j-1;
//                                    elements[ix]=elements[ix]-w[i-1]*e[j-1]/e[kk];
//                                }
//                            }
//                        }
	                
//                        for (i=kk+1; i<=n; i++)
//                            mtxV.elements[(i-1)*n+kk-1]=e[i-1];
//                    }
//                }
//            }
	    
//            mm=n;
//            if (m+1<n) 
//                mm=m+1;
//            if (k<n) 
//                s[k]=elements[k*n+k];
//            if (m<mm) 
//                s[mm-1]=0.0;
//            if (l+1<mm) 
//                e[l]=elements[l*n+mm-1];

//            e[mm-1]=0.0;
//            nn=m;
//            if (m>n) 
//                nn=n;
//            if (nn>=k+1)
//            { 
//                for (j=k+1; j<=nn; j++)
//                { 
//                    for (i=1; i<=m; i++)
//                        mtxU.elements[(i-1)*m+j-1]=0.0;
//                    mtxU.elements[(j-1)*m+j-1]=1.0;
//                }
//            }
	    
//            if (k>=1)
//            { 
//                for (ll=1; ll<=k; ll++)
//                { 
//                    kk=k-ll+1; 
//                    iz=(kk-1)*m+kk-1;
//                    if (s[kk-1]!=0.0)
//                    { 
//                        if (nn>=kk+1)
//                        {
//                            for (j=kk+1; j<=nn; j++)
//                            { 
//                                d=0.0;
//                                for (i=kk; i<=m; i++)
//                                { 
//                                    ix=(i-1)*m+kk-1;
//                                    iy=(i-1)*m+j-1;
//                                    d=d+mtxU.elements[ix]*mtxU.elements[iy]/mtxU.elements[iz];
//                                }

//                                d=-d;
//                                for (i=kk; i<=m; i++)
//                                { 
//                                    ix=(i-1)*m+j-1;
//                                    iy=(i-1)*m+kk-1;
//                                    mtxU.elements[ix]=mtxU.elements[ix]+d*mtxU.elements[iy];
//                                }
//                            }
//                        }
	                  
//                        for (i=kk; i<=m; i++)
//                        { 
//                            ix=(i-1)*m+kk-1; 
//                            mtxU.elements[ix]=-mtxU.elements[ix];
//                        }

//                        mtxU.elements[iz]=1.0+mtxU.elements[iz];
//                        if (kk-1>=1)
//                        {
//                            for (i=1; i<=kk-1; i++)
//                                mtxU.elements[(i-1)*m+kk-1]=0.0;
//                        }
//                    }
//                    else
//                    { 
//                        for (i=1; i<=m; i++)
//                            mtxU.elements[(i-1)*m+kk-1]=0.0;
//                        mtxU.elements[(kk-1)*m+kk-1]=1.0;
//                    }
//                }
//            }

//            for (ll=1; ll<=n; ll++)
//            { 
//                kk=n-ll+1; 
//                iz=kk*n+kk-1;
	        
//                if ((kk<=l) && (e[kk-1]!=0.0))
//                { 
//                    for (j=kk+1; j<=n; j++)
//                    { 
//                        d=0.0;
//                        for (i=kk+1; i<=n; i++)
//                        { 
//                            ix=(i-1)*n+kk-1; 
//                            iy=(i-1)*n+j-1;
//                            d=d+mtxV.elements[ix]*mtxV.elements[iy]/mtxV.elements[iz];
//                        }
	                
//                        d=-d;
//                        for (i=kk+1; i<=n; i++)
//                        { 
//                            ix=(i-1)*n+j-1; 
//                            iy=(i-1)*n+kk-1;
//                            mtxV.elements[ix]=mtxV.elements[ix]+d*mtxV.elements[iy];
//                        }
//                    }
//                }
	        
//                for (i=1; i<=n; i++)
//                    mtxV.elements[(i-1)*n+kk-1]=0.0;
	        
//                mtxV.elements[iz-n]=1.0;
//            }
	    
//            for (i=1; i<=m; i++)
//                for (j=1; j<=n; j++)
//                    elements[(i-1)*n+j-1]=0.0;
	    
//            m1=mm; 
//            it=60;
//            while (true)
//            { 
//                if (mm==0)
//                { 
//                    ppp(elements,e,s,mtxV.elements,m,n);
//                    return true;
//                }
//                if (it==0)
//                { 
//                    ppp(elements,e,s,mtxV.elements,m,n);
//                    return false;
//                }
	        
//                kk=mm-1;
//                while ((kk!=0) && (Math.Abs(e[kk-1])!=0.0))
//                { 
//                    d=Math.Abs(s[kk-1])+Math.Abs(s[kk]);
//                    dd=Math.Abs(e[kk-1]);
//                    if (dd>eps*d) 
//                        kk=kk-1;
//                    else 
//                        e[kk-1]=0.0;
//                }
	        
//                if (kk==mm-1)
//                { 
//                    kk=kk+1;
//                    if (s[kk-1]<0.0)
//                    { 
//                        s[kk-1]=-s[kk-1];
//                        for (i=1; i<=n; i++)
//                        { 
//                            ix=(i-1)*n+kk-1; 
//                            mtxV.elements[ix]=-mtxV.elements[ix];}
//                    }
					
//                    while ((kk!=m1) && (s[kk-1]<s[kk]))
//                    { 
//                        d=s[kk-1]; 
//                        s[kk-1]=s[kk]; 
//                        s[kk]=d;
//                        if (kk<n)
//                        {
//                            for (i=1; i<=n; i++)
//                            { 
//                                ix=(i-1)*n+kk-1; 
//                                iy=(i-1)*n+kk;
//                                d=mtxV.elements[ix]; 
//                                mtxV.elements[ix]=mtxV.elements[iy]; 
//                                mtxV.elements[iy]=d;
//                            }
//                        }

//                        if (kk<m)
//                        {
//                            for (i=1; i<=m; i++)
//                            { 
//                                ix=(i-1)*m+kk-1; 
//                                iy=(i-1)*m+kk;
//                                d=mtxU.elements[ix]; 
//                                mtxU.elements[ix]=mtxU.elements[iy]; 
//                                mtxU.elements[iy]=d;
//                            }
//                        }

//                        kk=kk+1;
//                    }
	            
//                    it=60;
//                    mm=mm-1;
//                }
//                else
//                { 
//                    ks=mm;
//                    while ((ks>kk) && (Math.Abs(s[ks-1])!=0.0))
//                    { 
//                        d=0.0;
//                        if (ks!=mm) 
//                            d=d+Math.Abs(e[ks-1]);
//                        if (ks!=kk+1) 
//                            d=d+Math.Abs(e[ks-2]);
	                
//                        dd=Math.Abs(s[ks-1]);
//                        if (dd>eps*d) 
//                            ks=ks-1;
//                        else 
//                            s[ks-1]=0.0;
//                    }
	            
//                    if (ks==kk)
//                    { 
//                        kk=kk+1;
//                        d=Math.Abs(s[mm-1]);
//                        t=Math.Abs(s[mm-2]);
//                        if (t>d) 
//                            d=t;
	                
//                        t=Math.Abs(e[mm-2]);
//                        if (t>d) 
//                            d=t;
	                
//                        t=Math.Abs(s[kk-1]);
//                        if (t>d) 
//                            d=t;
	                
//                        t=Math.Abs(e[kk-1]);
//                        if (t>d) 
//                            d=t;
	                
//                        sm=s[mm-1]/d; 
//                        sm1=s[mm-2]/d;
//                        em1=e[mm-2]/d;
//                        sk=s[kk-1]/d; 
//                        ek=e[kk-1]/d;
//                        b=((sm1+sm)*(sm1-sm)+em1*em1)/2.0;
//                        c=sm*em1; 
//                        c=c*c; 
//                        shh=0.0;

//                        if ((b!=0.0)||(c!=0.0))
//                        { 
//                            shh=Math.Sqrt(b*b+c);
//                            if (b<0.0) 
//                                shh=-shh;

//                            shh=c/(b+shh);
//                        }
	                
//                        fg[0]=(sk+sm)*(sk-sm)-shh;
//                        fg[1]=sk*ek;
//                        for (i=kk; i<=mm-1; i++)
//                        { 
//                            sss(fg,cs);
//                            if (i!=kk) 
//                                e[i-2]=fg[0];

//                            fg[0]=cs[0]*s[i-1]+cs[1]*e[i-1];
//                            e[i-1]=cs[0]*e[i-1]-cs[1]*s[i-1];
//                            fg[1]=cs[1]*s[i];
//                            s[i]=cs[0]*s[i];

//                            if ((cs[0]!=1.0)||(cs[1]!=0.0))
//                            {
//                                for (j=1; j<=n; j++)
//                                { 
//                                    ix=(j-1)*n+i-1;
//                                    iy=(j-1)*n+i;
//                                    d=cs[0]*mtxV.elements[ix]+cs[1]*mtxV.elements[iy];
//                                    mtxV.elements[iy]=-cs[1]*mtxV.elements[ix]+cs[0]*mtxV.elements[iy];
//                                    mtxV.elements[ix]=d;
//                                }
//                            }

//                            sss(fg,cs);
//                            s[i-1]=fg[0];
//                            fg[0]=cs[0]*e[i-1]+cs[1]*s[i];
//                            s[i]=-cs[1]*e[i-1]+cs[0]*s[i];
//                            fg[1]=cs[1]*e[i];
//                            e[i]=cs[0]*e[i];

//                            if (i<m)
//                            {
//                                if ((cs[0]!=1.0)||(cs[1]!=0.0))
//                                {
//                                    for (j=1; j<=m; j++)
//                                    { 
//                                        ix=(j-1)*m+i-1;
//                                        iy=(j-1)*m+i;
//                                        d=cs[0]*mtxU.elements[ix]+cs[1]*mtxU.elements[iy];
//                                        mtxU.elements[iy]=-cs[1]*mtxU.elements[ix]+cs[0]*mtxU.elements[iy];
//                                        mtxU.elements[ix]=d;
//                                    }
//                                }
//                            }
//                        }
	                
//                        e[mm-2]=fg[0];
//                        it=it-1;
//                    }
//                    else
//                    { 
//                        if (ks==mm)
//                        { 
//                            kk=kk+1;
//                            fg[1]=e[mm-2]; 
//                            e[mm-2]=0.0;
//                            for (ll=kk; ll<=mm-1; ll++)
//                            { 
//                                i=mm+kk-ll-1;
//                                fg[0]=s[i-1];
//                                sss(fg,cs);
//                                s[i-1]=fg[0];
//                                if (i!=kk)
//                                { 
//                                    fg[1]=-cs[1]*e[i-2];
//                                    e[i-2]=cs[0]*e[i-2];
//                                }
	                        
//                                if ((cs[0]!=1.0)||(cs[1]!=0.0))
//                                {
//                                    for (j=1; j<=n; j++)
//                                    { 
//                                        ix=(j-1)*n+i-1;
//                                        iy=(j-1)*n+mm-1;
//                                        d=cs[0]*mtxV.elements[ix]+cs[1]*mtxV.elements[iy];
//                                        mtxV.elements[iy]=-cs[1]*mtxV.elements[ix]+cs[0]*mtxV.elements[iy];
//                                        mtxV.elements[ix]=d;
//                                    }
//                                }
//                            }
//                        }
//                        else
//                        { 
//                            kk=ks+1;
//                            fg[1]=e[kk-2];
//                            e[kk-2]=0.0;
//                            for (i=kk; i<=mm; i++)
//                            { 
//                                fg[0]=s[i-1];
//                                sss(fg,cs);
//                                s[i-1]=fg[0];
//                                fg[1]=-cs[1]*e[i-1];
//                                e[i-1]=cs[0]*e[i-1];
//                                if ((cs[0]!=1.0)||(cs[1]!=0.0))
//                                {
//                                    for (j=1; j<=m; j++)
//                                    { 
//                                        ix=(j-1)*m+i-1;
//                                        iy=(j-1)*m+kk-2;
//                                        d=cs[0]*mtxU.elements[ix]+cs[1]*mtxU.elements[iy];
//                                        mtxU.elements[iy]=-cs[1]*mtxU.elements[ix]+cs[0]*mtxU.elements[iy];
//                                        mtxU.elements[ix]=d;
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }

//        /**
//         * 内部函数，由SplitUV函数调用
//         */
//        private void ppp(double[] a, double[] e, double[] s, double[] v, int m, int n)
//        { 
//            int i,j,p,q;
//            double d;

//            if (m>=n) 
//                i=n;
//            else 
//                i=m;

//            for (j=1; j<=i-1; j++)
//            { 
//                a[(j-1)*n+j-1]=s[j-1];
//                a[(j-1)*n+j]=e[j-1];
//            }
	    
//            a[(i-1)*n+i-1]=s[i-1];
//            if (m<n) 
//                a[(i-1)*n+i]=e[i-1];
	    
//            for (i=1; i<=n-1; i++)
//            {
//                for (j=i+1; j<=n; j++)
//                { 
//                    p=(i-1)*n+j-1; 
//                    q=(j-1)*n+i-1;
//                    d=v[p]; 
//                    v[p]=v[q]; 
//                    v[q]=d;
//                }
//            }
//        }

//        /**
//         * 内部函数，由SplitUV函数调用
//         */
//        private void sss(double[] fg, double[] cs)
//        { 
//            double r,d;
	    
//            if ((Math.Abs(fg[0])+Math.Abs(fg[1]))==0.0)
//            { 
//                cs[0]=1.0; 
//                cs[1]=0.0; 
//                d=0.0;
//            }
//            else 
//            { 
//                d=Math.Sqrt(fg[0]*fg[0]+fg[1]*fg[1]);
//                if (Math.Abs(fg[0])>Math.Abs(fg[1]))
//                { 
//                    d=Math.Abs(d);
//                    if (fg[0]<0.0) 
//                        d=-d;
//                }
//                if (Math.Abs(fg[1])>=Math.Abs(fg[0]))
//                { 
//                    d=Math.Abs(d);
//                    if (fg[1]<0.0) 
//                        d=-d;
//                }
	        
//                cs[0]=fg[0]/d; 
//                cs[1]=fg[1]/d;
//            }
	    
//            r=1.0;
//            if (Math.Abs(fg[0])>Math.Abs(fg[1])) 
//                r=cs[1];
//            else if (cs[0]!=0.0) 
//                r=1.0/cs[0];

//            fg[0]=d; 
//            fg[1]=r;
//        }

//        /**
//         * 求广义逆的奇异值分解法，分解成功后，原矩阵对角线元素就是矩阵的奇异值
//         * 
//         * @param mtxAP - 返回原矩阵的广义逆矩阵
//         * @param mtxU - 返回分解后的U矩阵
//         * @param mtxV - 返回分解后的V矩阵
//         * @param eps - 计算精度
//         * @return bool型，求解是否成功
//         */
//        public bool InvertUV(Matrix mtxAP, Matrix mtxU, Matrix mtxV, double eps)
//        { 
//            int i,j,k,l,t,p,q,f;

//            // 调用奇异值分解
//            if (! SplitUV(mtxU, mtxV, eps))
//                return false;

//            int m = numRows;
//            int n = numColumns;

//            // 初始化广义逆矩阵
//            if (! mtxAP.Init(n, m))
//                return false;

//            // 计算广义逆矩阵

//            j=n;
//            if (m<n) 
//                j=m;
//            j=j-1;
//            k=0;
//            while ((k<=j) && (elements[k*n+k]!=0.0)) 
//                k=k+1;

//            k=k-1;
//            for (i=0; i<=n-1; i++)
//            {
//                for (j=0; j<=m-1; j++)
//                { 
//                    t=i*m+j;	
//                    mtxAP.elements[t]=0.0;
//                    for (l=0; l<=k; l++)
//                    { 
//                        f=l*n+i; 
//                        p=j*m+l; 
//                        q=l*n+l;
//                        mtxAP.elements[t]=mtxAP.elements[t]+mtxV.elements[f]*mtxU.elements[p]/elements[q];
//                    }
//                }
//            }

//            return true;
//        }

//        /**
//         * 约化对称矩阵为对称三对角阵的豪斯荷尔德变换法
//         * 
//         * @param mtxQ - 返回豪斯荷尔德变换的乘积矩阵Q
//         * @param mtxT - 返回求得的对称三对角阵
//         * @param dblB - 一维数组，长度为矩阵的阶数，返回对称三对角阵的主对角线元素
//         * @param dblC - 一维数组，长度为矩阵的阶数，前n-1个元素返回对称三对角阵的
//         *               次对角线元素
//         * @return bool型，求解是否成功
//         */
//        public bool MakeSymTri(Matrix mtxQ, Matrix mtxT, double[] dblB, double[] dblC)
//        { 
//            int i,j,k,u;
//            double h,f,g,h2;
	    
//            // 初始化矩阵Q和T
//            if (! mtxQ.Init(numColumns, numColumns) ||
//                ! mtxT.Init(numColumns, numColumns))
//                return false;

//            if (dblB == null || dblC == null)
//                return false;

//            for (i=0; i<=numColumns-1; i++)
//            {
//                for (j=0; j<=numColumns-1; j++)
//                { 
//                    u=i*numColumns+j; 
//                    mtxQ.elements[u]=elements[u];
//                }
//            }

//            for (i=numColumns-1; i>=1; i--)
//            { 
//                h=0.0;
//                if (i>1)
//                {
//                    for (k=0; k<=i-1; k++)
//                    { 
//                        u=i*numColumns+k; 
//                        h=h+mtxQ.elements[u]*mtxQ.elements[u];
//                    }
//                }

//                if (h == 0.0)
//                { 
//                    dblC[i]=0.0;
//                    if (i==1) 
//                        dblC[i]=mtxQ.elements[i*numColumns+i-1];
//                    dblB[i]=0.0;
//                }
//                else
//                { 
//                    dblC[i]=Math.Sqrt(h);
//                    u=i*numColumns+i-1;
//                    if (mtxQ.elements[u]>0.0) 
//                        dblC[i]=-dblC[i];

//                    h=h-mtxQ.elements[u]*dblC[i];
//                    mtxQ.elements[u]=mtxQ.elements[u]-dblC[i];
//                    f=0.0;
//                    for (j=0; j<=i-1; j++)
//                    { 
//                        mtxQ.elements[j*numColumns+i]=mtxQ.elements[i*numColumns+j]/h;
//                        g=0.0;
//                        for (k=0; k<=j; k++)
//                            g=g+mtxQ.elements[j*numColumns+k]*mtxQ.elements[i*numColumns+k];

//                        if (j+1<=i-1)
//                            for (k=j+1; k<=i-1; k++)
//                                g=g+mtxQ.elements[k*numColumns+j]*mtxQ.elements[i*numColumns+k];

//                        dblC[j]=g/h;
//                        f=f+g*mtxQ.elements[j*numColumns+i];
//                    }
	            
//                    h2=f/(h+h);
//                    for (j=0; j<=i-1; j++)
//                    { 
//                        f=mtxQ.elements[i*numColumns+j];
//                        g=dblC[j]-h2*f;
//                        dblC[j]=g;
//                        for (k=0; k<=j; k++)
//                        { 
//                            u=j*numColumns+k;
//                            mtxQ.elements[u]=mtxQ.elements[u]-f*dblC[k]-g*mtxQ.elements[i*numColumns+k];
//                        }
//                    }
	            
//                    dblB[i]=h;
//                }
//            }
	    
//            for (i=0; i<=numColumns-2; i++) 
//                dblC[i]=dblC[i+1];
	    
//            dblC[numColumns-1]=0.0;
//            dblB[0]=0.0;
//            for (i=0; i<=numColumns-1; i++)
//            { 
//                if ((dblB[i]!=(double)0.0) && (i-1>=0))
//                {
//                    for (j=0; j<=i-1; j++)
//                    { 
//                        g=0.0;
//                        for (k=0; k<=i-1; k++)
//                            g=g+mtxQ.elements[i*numColumns+k]*mtxQ.elements[k*numColumns+j];

//                        for (k=0; k<=i-1; k++)
//                        { 
//                            u=k*numColumns+j;
//                            mtxQ.elements[u]=mtxQ.elements[u]-g*mtxQ.elements[k*numColumns+i];
//                        }
//                    }
//                }

//                u=i*numColumns+i;
//                dblB[i]=mtxQ.elements[u]; mtxQ.elements[u]=1.0;
//                if (i-1>=0)
//                {
//                    for (j=0; j<=i-1; j++)
//                    { 
//                        mtxQ.elements[i*numColumns+j]=0.0; 
//                        mtxQ.elements[j*numColumns+i]=0.0;
//                    }
//                }
//            }

//            // 构造对称三对角矩阵
//            for (i=0; i<numColumns; ++i)
//            {
//                for (j=0; j<numColumns; ++j)
//                {
//                    mtxT.SetElement(i, j, 0);
//                    k = i - j;
//                    if (k == 0) 
//                        mtxT.SetElement(i, j, dblB[j]);
//                    else if (k == 1)
//                        mtxT.SetElement(i, j, dblC[j]);
//                    else if (k == -1)
//                        mtxT.SetElement(i, j, dblC[i]);
//                }
//            }

//            return true;
//        }

//        /**
//         * 实对称三对角阵的全部特征值与特征向量的计算
//         * 
//         * @param dblB - 一维数组，长度为矩阵的阶数，传入对称三对角阵的主对角线元素；
//         *			     返回时存放全部特征值。
//         * @param dblC - 一维数组，长度为矩阵的阶数，前n-1个元素传入对称三对角阵的
//         *               次对角线元素
//         * @param mtxQ - 如果传入单位矩阵，则返回实对称三对角阵的特征值向量矩阵；
//         *			     如果传入MakeSymTri函数求得的矩阵A的豪斯荷尔德变换的乘积
//         *               矩阵Q，则返回矩阵A的特征值向量矩阵。其中第i列为与数组dblB
//         *               中第j个特征值对应的特征向量。
//         * @param nMaxIt - 迭代次数
//         * @param eps - 计算精度
//         * @return bool型，求解是否成功
//         */
//        public bool ComputeEvSymTri(double[] dblB, double[] dblC, Matrix mtxQ, int nMaxIt, double eps)
//        {
//            int i,j,k,m,it,u,v;
//            double d,f,h,g,p,r,e,s;
	    
//            // 初值
//            int n = mtxQ.GetNumColumns();
//            dblC[n-1]=0.0; 
//            d=0.0; 
//            f=0.0;
	    
//            // 迭代计算

//            for (j=0; j<=n-1; j++)
//            { 
//                it=0;
//                h=eps*(Math.Abs(dblB[j])+Math.Abs(dblC[j]));
//                if (h>d) 
//                    d=h;
	        
//                m=j;
//                while ((m<=n-1) && (Math.Abs(dblC[m])>d)) 
//                    m=m+1;
	        
//                if (m!=j)
//                { 
//                    do
//                    { 
//                        if (it==nMaxIt)
//                            return false;

//                        it=it+1;
//                        g=dblB[j];
//                        p=(dblB[j+1]-g)/(2.0*dblC[j]);
//                        r=Math.Sqrt(p*p+1.0);
//                        if (p>=0.0) 
//                            dblB[j]=dblC[j]/(p+r);
//                        else 
//                            dblB[j]=dblC[j]/(p-r);
	                
//                        h=g-dblB[j];
//                        for (i=j+1; i<=n-1; i++)
//                            dblB[i]=dblB[i]-h;
	                
//                        f=f+h; 
//                        p=dblB[m]; 
//                        e=1.0; 
//                        s=0.0;
//                        for (i=m-1; i>=j; i--)
//                        { 
//                            g=e*dblC[i]; 
//                            h=e*p;
//                            if (Math.Abs(p)>=Math.Abs(dblC[i]))
//                            { 
//                                e=dblC[i]/p; 
//                                r=Math.Sqrt(e*e+1.0);
//                                dblC[i+1]=s*p*r; 
//                                s=e/r; 
//                                e=1.0/r;
//                            }
//                            else
//                            { 
//                                e=p/dblC[i]; 
//                                r=Math.Sqrt(e*e+1.0);
//                                dblC[i+1]=s*dblC[i]*r;
//                                s=1.0/r; 
//                                e=e/r;
//                            }
	                    
//                            p=e*dblB[i]-s*g;
//                            dblB[i+1]=h+s*(e*g+s*dblB[i]);
//                            for (k=0; k<=n-1; k++)
//                            { 
//                                u=k*n+i+1; 
//                                v=u-1;
//                                h=mtxQ.elements[u]; 
//                                mtxQ.elements[u]=s*mtxQ.elements[v]+e*h;
//                                mtxQ.elements[v]=e*mtxQ.elements[v]-s*h;
//                            }
//                        }
	                
//                        dblC[j]=s*p; 
//                        dblB[j]=e*p;
	            
//                    } while (Math.Abs(dblC[j])>d);
//                }
	        
//                dblB[j]=dblB[j]+f;
//            }
	    
//            for (i=0; i<=n-1; i++)
//            { 
//                k=i; 
//                p=dblB[i];
//                if (i+1<=n-1)
//                { 
//                    j=i+1;
//                    while ((j<=n-1) && (dblB[j]<=p))
//                    { 
//                        k=j; 
//                        p=dblB[j]; 
//                        j=j+1;
//                    }
//                }

//                if (k!=i)
//                { 
//                    dblB[k]=dblB[i]; 
//                    dblB[i]=p;
//                    for (j=0; j<=n-1; j++)
//                    { 
//                        u=j*n+i; 
//                        v=j*n+k;
//                        p=mtxQ.elements[u]; 
//                        mtxQ.elements[u]=mtxQ.elements[v]; 
//                        mtxQ.elements[v]=p;
//                    }
//                }
//            }
	    
//            return true;
//        }

//        /**
//         * 约化一般实矩阵为赫申伯格矩阵的初等相似变换法
//         */
//        public void MakeHberg()
//        { 
//            int i = 0,j,k,u,v;
//            double d,t;

//            for (k=1; k<=numColumns-2; k++)
//            { 
//                d=0.0;
//                for (j=k; j<=numColumns-1; j++)
//                { 
//                    u=j*numColumns+k-1; 
//                    t=elements[u];
//                    if (Math.Abs(t)>Math.Abs(d))
//                    { 
//                        d=t; 
//                        i=j;
//                    }
//                }
	        
//                if (d != 0.0)
//                { 
//                    if (i!=k)
//                    { 
//                        for (j=k-1; j<=numColumns-1; j++)
//                        { 
//                            u=i*numColumns+j; 
//                            v=k*numColumns+j;
//                            t=elements[u]; 
//                            elements[u]=elements[v]; 
//                            elements[v]=t;
//                        }
	                
//                        for (j=0; j<=numColumns-1; j++)
//                        { 
//                            u=j*numColumns+i; 
//                            v=j*numColumns+k;
//                            t=elements[u]; 
//                            elements[u]=elements[v]; 
//                            elements[v]=t;
//                        }
//                    }
	            
//                    for (i=k+1; i<=numColumns-1; i++)
//                    { 
//                        u=i*numColumns+k-1; 
//                        t=elements[u]/d; 
//                        elements[u]=0.0;
//                        for (j=k; j<=numColumns-1; j++)
//                        { 
//                            v=i*numColumns+j;
//                            elements[v]=elements[v]-t*elements[k*numColumns+j];
//                        }
	                
//                        for (j=0; j<=numColumns-1; j++)
//                        { 
//                            v=j*numColumns+k;
//                            elements[v]=elements[v]+t*elements[j*numColumns+i];
//                        }
//                    }
//                }
//            }
//        }

//        /**
//         * 求赫申伯格矩阵全部特征值的QR方法
//         * 
//         * @param dblU - 一维数组，长度为矩阵的阶数，返回时存放特征值的实部
//         * @param dblV - 一维数组，长度为矩阵的阶数，返回时存放特征值的虚部
//         * @param nMaxIt - 迭代次数
//         * @param eps - 计算精度
//         * @return bool型，求解是否成功
//         */
//        public bool ComputeEvHBerg(double[] dblU, double[] dblV, int nMaxIt, double eps)
//        { 
//            int m,it,i,j,k,l,ii,jj,kk,ll;
//            double b,c,w,g,xy,p,q,r,x,s,e,f,z,y;
	    
//            int n = numColumns;

//            it=0; 
//            m=n;
//            while (m!=0)
//            { 
//                l=m-1;
//                while ((l>0) && (Math.Abs(elements[l*n+l-1]) > 
//                    eps*(Math.Abs(elements[(l-1)*n+l-1])+Math.Abs(elements[l*n+l])))) 
//                    l=l-1;

//                ii=(m-1)*n+m-1; 
//                jj=(m-1)*n+m-2;
//                kk=(m-2)*n+m-1; 
//                ll=(m-2)*n+m-2;
//                if (l==m-1)
//                { 
//                    dblU[m-1]=elements[(m-1)*n+m-1]; 
//                    dblV[m-1]=0.0;
//                    m=m-1; 
//                    it=0;
//                }
//                else if (l==m-2)
//                { 
//                    b=-(elements[ii]+elements[ll]);
//                    c=elements[ii]*elements[ll]-elements[jj]*elements[kk];
//                    w=b*b-4.0*c;
//                    y=Math.Sqrt(Math.Abs(w));
//                    if (w>0.0)
//                    { 
//                        xy=1.0;
//                        if (b<0.0) 
//                            xy=-1.0;
//                        dblU[m-1]=(-b-xy*y)/2.0;
//                        dblU[m-2]=c/dblU[m-1];
//                        dblV[m-1]=0.0; dblV[m-2]=0.0;
//                    }
//                    else
//                    { 
//                        dblU[m-1]=-b/2.0; 
//                        dblU[m-2]=dblU[m-1];
//                        dblV[m-1]=y/2.0; 
//                        dblV[m-2]=-dblV[m-1];
//                    }
	            
//                    m=m-2; 
//                    it=0;
//                }
//                else
//                { 
//                    if (it>=nMaxIt)
//                        return false;

//                    it=it+1;
//                    for (j=l+2; j<=m-1; j++)
//                        elements[j*n+j-2]=0.0;
//                    for (j=l+3; j<=m-1; j++)
//                        elements[j*n+j-3]=0.0;
//                    for (k=l; k<=m-2; k++)
//                    { 
//                        if (k!=l)
//                        { 
//                            p=elements[k*n+k-1]; 
//                            q=elements[(k+1)*n+k-1];
//                            r=0.0;
//                            if (k!=m-2) 
//                                r=elements[(k+2)*n+k-1];
//                        }
//                        else
//                        { 
//                            x=elements[ii]+elements[ll];
//                            y=elements[ll]*elements[ii]-elements[kk]*elements[jj];
//                            ii=l*n+l; 
//                            jj=l*n+l+1;
//                            kk=(l+1)*n+l; 
//                            ll=(l+1)*n+l+1;
//                            p=elements[ii]*(elements[ii]-x)+elements[jj]*elements[kk]+y;
//                            q=elements[kk]*(elements[ii]+elements[ll]-x);
//                            r=elements[kk]*elements[(l+2)*n+l+1];
//                        }
	                
//                        if ((Math.Abs(p)+Math.Abs(q)+Math.Abs(r))!=0.0)
//                        { 
//                            xy=1.0;
//                            if (p<0.0) 
//                                xy=-1.0;
//                            s=xy*Math.Sqrt(p*p+q*q+r*r);
//                            if (k!=l) 
//                                elements[k*n+k-1]=-s;
//                            e=-q/s; 
//                            f=-r/s; 
//                            x=-p/s;
//                            y=-x-f*r/(p+s);
//                            g=e*r/(p+s);
//                            z=-x-e*q/(p+s);
//                            for (j=k; j<=m-1; j++)
//                            { 
//                                ii=k*n+j; 
//                                jj=(k+1)*n+j;
//                                p=x*elements[ii]+e*elements[jj];
//                                q=e*elements[ii]+y*elements[jj];
//                                r=f*elements[ii]+g*elements[jj];
//                                if (k!=m-2)
//                                { 
//                                    kk=(k+2)*n+j;
//                                    p=p+f*elements[kk];
//                                    q=q+g*elements[kk];
//                                    r=r+z*elements[kk]; 
//                                    elements[kk]=r;
//                                }
	                        
//                                elements[jj]=q; elements[ii]=p;
//                            }
	                    
//                            j=k+3;
//                            if (j>=m-1) 
//                                j=m-1;
	                    
//                            for (i=l; i<=j; i++)
//                            { 
//                                ii=i*n+k; 
//                                jj=i*n+k+1;
//                                p=x*elements[ii]+e*elements[jj];
//                                q=e*elements[ii]+y*elements[jj];
//                                r=f*elements[ii]+g*elements[jj];
//                                if (k!=m-2)
//                                { 
//                                    kk=i*n+k+2;
//                                    p=p+f*elements[kk];
//                                    q=q+g*elements[kk];
//                                    r=r+z*elements[kk]; 
//                                    elements[kk]=r;
//                                }
	                        
//                                elements[jj]=q; 
//                                elements[ii]=p;
//                            }
//                        }
//                    }
//                }
//            }
	    
//            return true;
//        }

//        /**
//         * 求实对称矩阵特征值与特征向量的雅可比法
//         * 
//         * @param dblEigenValue - 一维数组，长度为矩阵的阶数，返回时存放特征值
//         * @param mtxEigenVector - 返回时存放特征向量矩阵，其中第i列为与数组
//         *                         dblEigenValue中第j个特征值对应的特征向量
//         * @param nMaxIt - 迭代次数
//         * @param eps - 计算精度
//         * @return bool型，求解是否成功
//         */
//        public bool ComputeEvJacobi(double[] dblEigenValue, Matrix mtxEigenVector, int nMaxIt, double eps)
//        { 
//            int i,j,p = 0,q = 0,u,w,t,s,l;
//            double fm,cn,sn,omega,x,y,d;
	    
//            if (! mtxEigenVector.Init(numColumns, numColumns))
//                return false;

//            l=1;
//            for (i=0; i<=numColumns-1; i++)
//            { 
//                mtxEigenVector.elements[i*numColumns+i]=1.0;
//                for (j=0; j<=numColumns-1; j++)
//                    if (i!=j) 
//                        mtxEigenVector.elements[i*numColumns+j]=0.0;
//            }
	    
//            while (true)
//            { 
//                fm=0.0;
//                for (i=1; i<=numColumns-1; i++)
//                {
//                    for (j=0; j<=i-1; j++)
//                    { 
//                        d=Math.Abs(elements[i*numColumns+j]);
//                        if ((i!=j) && (d>fm))
//                        { 
//                            fm=d; 
//                            p=i; 
//                            q=j;
//                        }
//                    }
//                }

//                if (fm<eps)
//                {
//                    for (i=0; i<numColumns; ++i)
//                        dblEigenValue[i] = GetElement(i,i);
//                    return true;
//                }

//                if (l>nMaxIt)  
//                    return false;
	        
//                l=l+1;
//                u=p*numColumns+q; 
//                w=p*numColumns+p; 
//                t=q*numColumns+p; 
//                s=q*numColumns+q;
//                x=-elements[u]; 
//                y=(elements[s]-elements[w])/2.0;
//                omega=x/Math.Sqrt(x*x+y*y);

//                if (y<0.0) 
//                    omega=-omega;

//                sn=1.0+Math.Sqrt(1.0-omega*omega);
//                sn=omega/Math.Sqrt(2.0*sn);
//                cn=Math.Sqrt(1.0-sn*sn);
//                fm=elements[w];
//                elements[w]=fm*cn*cn+elements[s]*sn*sn+elements[u]*omega;
//                elements[s]=fm*sn*sn+elements[s]*cn*cn-elements[u]*omega;
//                elements[u]=0.0; 
//                elements[t]=0.0;
//                for (j=0; j<=numColumns-1; j++)
//                {
//                    if ((j!=p) && (j!=q))
//                    { 
//                        u=p*numColumns+j; w=q*numColumns+j;
//                        fm=elements[u];
//                        elements[u]=fm*cn+elements[w]*sn;
//                        elements[w]=-fm*sn+elements[w]*cn;
//                    }
//                }

//                for (i=0; i<=numColumns-1; i++)
//                {
//                    if ((i!=p) && (i!=q))
//                    { 
//                        u=i*numColumns+p; 
//                        w=i*numColumns+q;
//                        fm=elements[u];
//                        elements[u]=fm*cn+elements[w]*sn;
//                        elements[w]=-fm*sn+elements[w]*cn;
//                    }
//                }

//                for (i=0; i<=numColumns-1; i++)
//                { 
//                    u=i*numColumns+p; 
//                    w=i*numColumns+q;
//                    fm=mtxEigenVector.elements[u];
//                    mtxEigenVector.elements[u]=fm*cn+mtxEigenVector.elements[w]*sn;
//                    mtxEigenVector.elements[w]=-fm*sn+mtxEigenVector.elements[w]*cn;
//                }
//            }
//        }

//        /**
//         * 求实对称矩阵特征值与特征向量的雅可比过关法
//         * 
//         * @param dblEigenValue - 一维数组，长度为矩阵的阶数，返回时存放特征值
//         * @param mtxEigenVector - 返回时存放特征向量矩阵，其中第i列为与数组
//         *                         dblEigenValue中第j个特征值对应的特征向量
//         * @param eps - 计算精度
//         * @return bool型，求解是否成功
//         */
//        public bool ComputeEvJacobi(double[] dblEigenValue, Matrix mtxEigenVector, double eps)
//        { 
//            int i,j,p,q,u,w,t,s;
//            double ff,fm,cn,sn,omega,x,y,d;
	    
//            if (! mtxEigenVector.Init(numColumns, numColumns))
//                return false;

//            for (i=0; i<=numColumns-1; i++)
//            { 
//                mtxEigenVector.elements[i*numColumns+i]=1.0;
//                for (j=0; j<=numColumns-1; j++)
//                    if (i!=j) 
//                        mtxEigenVector.elements[i*numColumns+j]=0.0;
//            }
	    
//            ff=0.0;
//            for (i=1; i<=numColumns-1; i++)
//            {
//                for (j=0; j<=i-1; j++)
//                { 
//                    d=elements[i*numColumns+j]; 
//                    ff=ff+d*d; 
//                }
//            }

//            ff=Math.Sqrt(2.0*ff);
//            ff=ff/(1.0*numColumns);

//            bool nextLoop = false;
//            while (true)
//            {
//                for (i=1; i<=numColumns-1; i++)
//                {
//                    for (j=0; j<=i-1; j++)
//                    { 
//                        d=Math.Abs(elements[i*numColumns+j]);
//                        if (d>ff)
//                        { 
//                            p=i; 
//                            q=j;

//                            u=p*numColumns+q; 
//                            w=p*numColumns+p; 
//                            t=q*numColumns+p; 
//                            s=q*numColumns+q;
//                            x=-elements[u]; 
//                            y=(elements[s]-elements[w])/2.0;
//                            omega=x/Math.Sqrt(x*x+y*y);
//                            if (y<0.0) 
//                                omega=-omega;
					    
//                            sn=1.0+Math.Sqrt(1.0-omega*omega);
//                            sn=omega/Math.Sqrt(2.0*sn);
//                            cn=Math.Sqrt(1.0-sn*sn);
//                            fm=elements[w];
//                            elements[w]=fm*cn*cn+elements[s]*sn*sn+elements[u]*omega;
//                            elements[s]=fm*sn*sn+elements[s]*cn*cn-elements[u]*omega;
//                            elements[u]=0.0; elements[t]=0.0;
					    
//                            for (j=0; j<=numColumns-1; j++)
//                            {
//                                if ((j!=p)&&(j!=q))
//                                { 
//                                    u=p*numColumns+j; 
//                                    w=q*numColumns+j;
//                                    fm=elements[u];
//                                    elements[u]=fm*cn+elements[w]*sn;
//                                    elements[w]=-fm*sn+elements[w]*cn;
//                                }
//                            }

//                            for (i=0; i<=numColumns-1; i++)
//                            {
//                                if ((i!=p)&&(i!=q))
//                                { 
//                                    u=i*numColumns+p; 
//                                    w=i*numColumns+q;
//                                    fm=elements[u];
//                                    elements[u]=fm*cn+elements[w]*sn;
//                                    elements[w]=-fm*sn+elements[w]*cn;
//                                }
//                            }
					    
//                            for (i=0; i<=numColumns-1; i++)
//                            { 
//                                u=i*numColumns+p; 
//                                w=i*numColumns+q;
//                                fm=mtxEigenVector.elements[u];
//                                mtxEigenVector.elements[u]=fm*cn+mtxEigenVector.elements[w]*sn;
//                                mtxEigenVector.elements[w]=-fm*sn+mtxEigenVector.elements[w]*cn;
//                            }

//                            nextLoop = true;
//                            break;
//                        }
//                    }

//                    if (nextLoop)
//                        break;
//                }
		        
//                if (nextLoop)
//                {
//                    nextLoop = false;
//                    continue;
//                }

//                nextLoop = false;

//                // 如果达到精度要求，退出循环，返回结果
//                if (ff<eps) 
//                {
//                    for (i=0; i<numColumns; ++i)
//                        dblEigenValue[i] = GetElement(i,i);
//                    return true;
//                }
		    
//                ff=ff/(1.0*numColumns);
//            }
//        }
//    }
//}
