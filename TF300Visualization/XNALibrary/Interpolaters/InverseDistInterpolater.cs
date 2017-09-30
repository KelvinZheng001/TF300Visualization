using System;
using System.Collections.Generic;
using System.Text;
using XNAHelper.Interpolaters;

namespace TF300.App.GUI.DatabaseUI.XNALibrary.Interpolaters
{
    public class InverseDistInterpolater : Interpolater
    {
        private int mExponent = 8;
        List<PointValue> mPoints;
        public InverseDistInterpolater(int exponent,List<PointValue> points)
        {
            mPoints = points;
            mExponent = exponent;
        }
        public double GetInterpolatedZ(float xpos, float ypos )
        {
            double z = 0;
            double totalW = 0;
            List<double> weightList = new List<double>();
            List<double> valueList = new List<double>();
            foreach (PointValue v in mPoints)
            {
                double distance = Math.Sqrt(Math.Pow(xpos - v.X, 2) + Math.Pow(ypos - v.Y, 2));
                if (Math.Abs(distance) > 0.00001)
                {
                    double weight = 1 / Math.Pow(distance, mExponent);
                    totalW += weight;
                    weightList.Add(weight);
                    valueList.Add(v.Value);
                }

            }
            for (int i = 0; i < weightList.Count; i++)
            {
                weightList[i] = weightList[i] / totalW;
            }
            for (int i = 0; i < weightList.Count; i++)
            {
                z += weightList[i] * valueList[i];
            }
            return z;
        }
    }
}
