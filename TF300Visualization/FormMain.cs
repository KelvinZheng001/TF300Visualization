using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Text;
using System.Windows.Forms;
using TF300.App.GUI.DatabaseUI.XNALibrary.Interpolaters;

namespace TF300Visualization
{
    public partial class FormMain : Form
    {
        public FormMain()
        {
            InitializeComponent();

            


            //this.Controls.Add
        }

        private XNAMesh3D mXna3D;

        private void button1_Click(object sender, EventArgs e)
        {
            mXna3D = new XNAMesh3D();

            mXna3D.Dock = DockStyle.Fill;
            this.GraphicPanel.Controls.Add(mXna3D);
            mXna3D.SendToBack();

            

            List<PointValue> inputData = new List<PointValue>();
            inputData.Add(new PointValue(10,10,0.02));
            inputData.Add(new PointValue(10,20,0.01));


            List<PointValue> _listProcessedMappingData = ProcessMappingData(inputData);
            mXna3D.SetData(_listProcessedMappingData, 150f);

        }


        Dictionary<PointF, List<PointValue>> _dictValueByPoint = new Dictionary<PointF, List<PointValue>>();
        /// <summary>
        /// 处理原始的点数据,处理重复点等信息.
        /// </summary>
        private List<PointValue> ProcessMappingData(List<PointValue> data)
        {

            _dictValueByPoint.Clear();
            //将所有点按点建索引字典
            foreach (PointValue pointvalue in data)
            {
                PointF newPoint = new PointF(pointvalue.X, pointvalue.Y);
                if (!_dictValueByPoint.ContainsKey(newPoint))
                {
                    _dictValueByPoint[newPoint] = new List<PointValue>();
                }
                _dictValueByPoint[newPoint].Add(pointvalue);
                //foreach (PointF point in _dictValueByPoint.Keys)
                //{
                //    if (pointvalue.X == point.X && pointvalue.Y == point.Y)
                //    {
                //        newPoint = point;
                //        _dictValueByPoint[point].Add(pointvalue.Value);
                //        break;
                //    }
                //}
                //if (newPoint == null)
                //{
                //    newPoint = new PointF(pointvalue.X, pointvalue.Y);
                //    _dictValueByPoint.Add(newPoint, new List<double>());
                //    _dictValueByPoint[newPoint].Add(pointvalue.Value);
                //}
            }
            List<PointValue> listResult = new List<PointValue>();
            //对每个点的数据求平均
            foreach (PointF point in _dictValueByPoint.Keys)
            {
                double sumValue = 0;
                foreach (PointValue pointValue in _dictValueByPoint[point])
                {
                    sumValue += pointValue.Value;
                }
                double averageValue = sumValue / _dictValueByPoint[point].Count;
                PointValue averagePointValue = new PointValue(point.X, point.Y, averageValue);
                listResult.Add(averagePointValue);
            }
            return listResult;
        }

        private void button2_Click(object sender, EventArgs e)
        {
            this.mXna3D.GrpahicType = GraphicTypeEnum.ContourLine3D;
        }

        private void button3_Click(object sender, EventArgs e)
        {
            this.mXna3D.GrpahicType = GraphicTypeEnum.ContourZone3D;
        }

        private void button4_Click(object sender, EventArgs e)
        {
            this.mXna3D.GrpahicType = GraphicTypeEnum.Color3D;
        }

        private void button5_Click(object sender, EventArgs e)
        {

        }
    }
}
