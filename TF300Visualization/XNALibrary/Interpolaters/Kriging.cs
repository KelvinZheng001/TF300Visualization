using System;
using System.Collections.Generic;
using System.Text;
using XNAHelper.Maths;
using TF300.App.GUI.DatabaseUI.XNALibrary.Interpolaters;

namespace XNAHelper.Interpolaters
{

    /// <summary>
    /// ����
    /// </summary>
    internal class Variogram
    {
        private double mBeta = 1.5;
        private double mAlpha;
        private double mNugsq;
        /// <summary>
        /// 
        /// </summary>
        /// <param name="points"></param>
        /// <param name="nug"></param>
        public Variogram(List<PointValue> points, double nug)
        {
            mNugsq = nug * nug;
            double rb = 0, nom = 0, denom = 0;
            for (int i = 0; i < points.Count; i++)
            {
                for (int j = i + 1; j < points.Count; j++)
                {
                    rb = Sqr(points[i], points[j]);
                    rb = Math.Pow(rb, 0.5 * mBeta);
                    nom += rb * (0.5 * (points[i].Value - points[j].Value) * (points[i].Value - points[j].Value) - mNugsq);
                    denom += rb * rb;
                }
            }
            mAlpha = nom / denom;
        }
        public double Vargram(double r)
        {
            return mAlpha * Math.Pow(r, mBeta);
        }
        /// <summary>
        /// ������x,y��ֵ��ƽ����
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        private double Sqr(PointValue a, PointValue b)
        {
            return Math.Pow(b.Y - a.Y, 2) + Math.Pow(b.X - a.X, 2);
        }

        private double GetDistance(PointValue a, PointValue b)
        {
            return Math.Sqrt(Math.Pow(b.Y - a.Y, 2) + Math.Pow(b.X - a.X, 2));
        }
    }
    public class Kriging : Interpolater
    {
        private int mSize = 0;
        private Matrix mAMatrix;
        private List<double> mlist = new List<double>();
        private double mSemivariance = 0;
        private Matrix mL;
        private Matrix mU;
        private Matrix mP;
        private Variogram mVariogram;

        private List<PointValue> mPoints = new List<PointValue>();
        /// <summary>
        /// ���������ֵ��
        /// </summary>
        /// <param name="points"></param>
        /// <param name="dSemivariance">�˲�����Ч</param>
        public Kriging(List<PointValue> points, double dSemivariance)
        {
            
            //����ʹ�õ�L/U�ⷽ�̷����ڽ⺬�жԳƵ��,������0��������п��ܳ�����ֵ�쳣,��˸����е����ݶ�����0.00001�Ա�����󲿷ֵ���ֵ������ֵĻ���
            //������Ӧ��ͨ��Ѱ�Ҹ��ȶ������Է�������ⷽ�������.

            mPoints = new List<PointValue>();
            Random random = new Random();
            foreach (PointValue point in points)
            {
                float randomValue = 0.05f * (float)random.NextDouble();
                PointValue value = new PointValue(point.X + randomValue, point.Y - randomValue, point.Value);
                mPoints.Add(value);
            }

            mSemivariance = dSemivariance;
            mSize = mPoints.Count;
            mVariogram = new Variogram(mPoints, 0.3);
            mAMatrix = new Matrix(mSize + 1, mSize + 1);

            mL = new Matrix(mSize + 1, mSize + 1);
            mU = new Matrix(mSize + 1, mSize + 1);
            mP = new Matrix(mSize + 1, mSize + 1);
            //��ʼ������
            for (int i = 1; i <= mSize + 1; i++)
            {
                for (int j = 1; j <= mSize + 1; j++)
                {
                    if (i == mSize + 1 && j == mSize + 1)
                    {
                        mAMatrix[i, j] = new Complex(0);
                    }
                    else if (i == mSize + 1 || j == mSize + 1)
                    {
                        mAMatrix[i, j] = new Complex(1);
                    }
                    else
                    {
                        double v = GetDistance(mPoints[i - 1], mPoints[j - 1]);//����
                        mAMatrix[i, j] = new Complex(mVariogram.Vargram(v));//����ķ���//new Complex(v * mSemivariance);
                    }
                }
            }

            GetParameters();//����P L U����  LU�ֽ� ׼��������Ȩ�ط������.
        }
        private void GetParameters()
        {
            Matrix A = mAMatrix.Clone();

            if (!A.IsSquare())
                throw new InvalidOperationException("Cannot uniquely solve non-square equation system.");

            mP = A.LUSafe();
            mL = A.ExtractLowerTrapeze() - Matrix.Diag(A.DiagVector()) + Matrix.Identity(A.RowCount);
            mU = A.ExtractUpperTrapeze();
        }

        private double GetDistance(PointValue a, PointValue b)
        {
            return Math.Sqrt(Math.Pow(b.Y - a.Y, 2) + Math.Pow(b.X - a.X, 2));
        }

        private double GetDistance(float x1, float y1, float x2, float y2)
        {
            return Math.Sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
        }

        /// <summary>
        /// ��ֵ
        /// </summary>
        /// <param name="xpos"></param>
        /// <param name="ypos"></param>
        /// <param name="points"></param>
        /// <returns></returns>
        public double GetInterpolatedZ(float xpos, float ypos)
        {
            Matrix b = new Matrix(mSize + 1, 1);
            for (int i = 1; i <= mSize; i++)
            {
                double dis = GetDistance(xpos, ypos, mPoints[i - 1].X, mPoints[i - 1].Y);
                b[i, 1] = new Complex(mVariogram.Vargram(dis));//new Complex(dis * mSemivariance);
            }
            b[mSize + 1, 1] = new Complex(1);

            b = mP * b;
            mL.ForwardInsertion(b);
            mU.BackwardInsertion(b);

            Matrix x = b;
            double Z = 0;
            double totalW1 = 0;
            for (int i = 1; i <= mSize; i++)
            {
                totalW1 += x[i, 1].Re;
            }
            for (int i = 1; i <= mSize; i++)
            {
                Z += mPoints[i - 1].Value * x[i, 1].Re / totalW1;
            }

            return Z;
        }
    }
}
