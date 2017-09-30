using System;
using System.Collections.Generic;
using System.Text;

namespace TF300.App.GUI.DatabaseUI.XNALibrary.Interpolaters
{
    public class PointValue
    {
        private float _X;

        public float X
        {
            get { return _X; }
            set { _X = value; }
        }
        private float _Y;

        public float Y
        {
            get { return _Y; }
            set { _Y = value; }
        }
        private double _Value;

        public double Value
        {
            get { return _Value; }
            set { _Value = value; }
        }
        public PointValue()
        {
        }
        public PointValue(float x, float y, double value)
        {
            this.X = x;
            this.Y = y;
            this.Value = value;
        }
    }

}
