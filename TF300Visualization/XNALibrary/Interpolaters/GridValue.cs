using System;
using System.Collections.Generic;
using System.Text;

namespace TF300.App.GUI.DatabaseUI.XNALibrary.Interpolaters
{
    public class GridValue
    {
        public int rows;
        public int collumns;
        public List<PointValue> ValueList;
        public GridValue()
        {
            rows = 0;
            collumns = 0;
            ValueList = new List<PointValue>();
        }
        public GridValue(int row, int collumn, List<PointValue> list)
        {
            rows = row;
            collumns = collumn;
            ValueList = list;
        }
    }

}
