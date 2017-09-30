using System;
using System.Collections.Generic;
using System.Text;
using Microsoft.Xna.Framework.Graphics;
using Microsoft.Xna.Framework;
using XNAHelper.Interpolaters;
using XNAHelper.CustomVertexs;
using XNAHelper.Helpers;
using TF300.App.GUI.DatabaseUI.XNALibrary.Interpolaters;
using System.Drawing.Imaging;
using System.Drawing;

namespace DrawUserPrimitives
{
    using Color = Microsoft.Xna.Framework.Graphics.Color;
    using ColorD = System.Drawing.Color;
    using System.Drawing.Drawing2D;
    //public class CustomMesh
    //{
    //    public VertexPositionNormalTexture[] arrVertex = null;
    //    public int[] arrIndex = null;

    //    /// <summary>
    //    /// 通过GridValue构造
    //    /// </summary>
    //    /// <param name="grid"></param>
    //    public CustomMesh(GridValue grid)
    //    {

    //        arrVertex = new VertexPositionNormalTexture[grid.collumns * grid.rows];
    //        for (int i = 0; i < grid.collumns; i++)
    //        {
    //            for (int j = 0; j < grid.rows; j++)
    //            {
    //                PointValue point = grid.ValueList[i * grid.collumns + j];
    //                arrVertex[i * grid.collumns + j] = new VertexPositionNormalTexture(new Vector3(point.X, point.Y, (float)point.Value), Vector3.Forward, Vector2.One);
    //            }
    //        }
    //        //生成三角形索引
    //        int width = grid.collumns;
    //        int height = grid.rows;
    //        List<int> listIndex = new List<int>();
    //        for (int i = 0; i < width - 1; i++)
    //        {
    //            for (int j = 0; j < height - 1; j++)
    //            {
    //                listIndex.Add(i * width + j);
    //                listIndex.Add(i * width + j + 1);
    //                listIndex.Add((i + 1) * width + j);
    //                listIndex.Add((i + 1) * width + j);
    //                listIndex.Add(i * width + j + 1);
    //                listIndex.Add((i + 1) * width + j + 1);
    //            }
    //        }
    //        arrIndex = listIndex.ToArray();
    //    }
    //}




    public class MeshGenerator
    {
        private static int InterpolateDimension = 20;
        public static int GridDimension = 50;
        /// <summary>
        /// 生成网格数据
        /// </summary>
        /// <returns></returns>
        private static GridValue InterpolateToGridData(List<PointValue> listPointValue, float radius)
        {
            Interpolater interpolater = null;
            if (listPointValue.Count > 1)
            {
                interpolater = new Kriging(listPointValue, 10);
            }
            //Interpolater interpolater = new InverseDistInterpolater(2);
            GridValue grid = new GridValue();
            grid.collumns = InterpolateDimension;
            grid.rows = InterpolateDimension;

            float start = -radius;
            float step = radius * 2 / (float)InterpolateDimension;
            //进行插值
            for (int i = 0; i < InterpolateDimension; i++)
            {
                for (int j = 0; j < InterpolateDimension; j++)
                {
                    PointValue pv = new PointValue();
                    if (listPointValue.Count > 1)
                    {
                        pv.Value = interpolater.GetInterpolatedZ(start + j * step, start + i * step);

                    }
                    else if (listPointValue.Count == 1)
                    {
                        pv.Value = listPointValue[0].Value;
                    }
                    pv.X = (float)j - InterpolateDimension / 2; pv.Y = (float)i - InterpolateDimension / 2;
                    grid.ValueList.Add(pv);

                }
            }
            return grid;
        }

        /// <summary>
        /// 转换数据点的坐标为3D图形坐标
        /// </summary>
        /// <param name="listPointValue"></param>
        /// <param name="radius"></param>
        /// <returns></returns>
        public static VertexPositionColorNormal[] TranslateDataPoints(List<PointValue> listPointValue, float radius)
        {
            VertexPositionColorNormal[] verticiesDataPoint = new VertexPositionColorNormal[listPointValue.Count];
            for (int i = 0; i < listPointValue.Count; i++)
            {
                float x = listPointValue[i].X * GridDimension / (radius * 2);
                float y = listPointValue[i].Y * GridDimension / (radius * 2);
                verticiesDataPoint[i] = new VertexPositionColorNormal(new Vector3(x, y, (float)listPointValue[i].Value), Color.Gray, Vector3.Up);
            }
            return verticiesDataPoint;

        }
        ///// <summary>
        ///// 生成网格数据,
        ///// </summary>
        ///// <param name="listPointValue">测量点数据</param>
        ///// <returns>返回适合XNA显示的数据格式</returns>
        //public static CustomMesh GetMesh(List<PointValue> listPointValue)
        //{
        //    GridValue grid = GetGridData(listPointValue);
        //    CustomMesh mesh = new CustomMesh(grid);
        //    return mesh;
        //}


        /// <summary>
        /// 用克里金法插值生成网格顶点
        /// </summary>
        /// <param name="listPointValue"></param>
        /// <param name="radius"></param>
        /// <returns></returns>
        public static VertexPositionColorNormal[] InterportationKrigingVertexPositionColor(List<PointValue> listPointValue, float radius)
        {


            GridValue grid = InterpolateToGridData(listPointValue, radius);//插值生成网格
            VertexPositionColorNormal[] vertices = new VertexPositionColorNormal[grid.collumns * grid.rows];


            for (int i = 0; i < grid.collumns; i++)
            {
                for (int j = 0; j < grid.rows; j++)
                {
                    PointValue point = grid.ValueList[i * grid.collumns + j];
                    vertices[i * grid.collumns + j] = new VertexPositionColorNormal(new Vector3(point.X, point.Y, (float)point.Value), Color.Red, Vector3.Zero);
                }
            }

            //增加网格密度
            //对原始散点数据先进行较低密度InterpolateDimension的克里金插值形成正交网格后,用CatmullRom插值法增加网络密度,以使得最终输出图像效果平滑
            VertexPositionColorNormal[] vertices1 = PopulateMeshCatmullRom(vertices, InterpolateDimension, GridDimension);
            return vertices1;


            //return vertices;
        }

        /// <summary>
        /// 网格增加加密度 使用CatmullRom插值
        /// </summary>
        /// <param name="sourceVerticies">源顶点</param>
        /// <param name="sourceDimension">源顶点的密度</param>
        /// <param name="destDimension">目标密度</param>
        /// <returns></returns>
        public static VertexPositionColorNormal[] PopulateMeshCatmullRom(VertexPositionColorNormal[] sourceVerticies, int sourceDimension, int destDimension)
        {
            VertexPositionColorNormal[] destVerticies = new VertexPositionColorNormal[destDimension * destDimension];

            //x,y是目标网格下标
            for (int y = 0; y < destDimension; y++)
            {
                for (int x = 0; x < destDimension; x++)
                {
                    //求源网格的取值点
                    int i = x * sourceDimension / destDimension;
                    int j = y * sourceDimension / destDimension;
                    //求插值点在x,y方向上的偏移量
                    float amountX = ((float)x / (float)destDimension - (float)i / (float)sourceDimension) / (1 / (float)sourceDimension);
                    float amountY = ((float)y / (float)destDimension - (float)j / (float)sourceDimension) / (1 / (float)sourceDimension);

                    float[] tempValue = new float[4];

                    float value0;
                    float value1;
                    float value2;
                    float value3;


                    for (int k = 0; k < 4; k++) //求偏移量为amountX的竖向四个临时插值点的值
                    {
                        int iMinus1 = i - 1 < 0 ? i : i - 1;
                        int iPlus1 = i + 1 > sourceDimension - 1 ? i : i + 1;
                        int iPlus2 = iPlus1 + 1 > sourceDimension - 1 ? iPlus1 : iPlus1 + 1;

                        int indexJ = j - 1 + k;
                        if (indexJ < 0)
                            indexJ = 0;
                        if (indexJ >= sourceDimension)
                            indexJ = sourceDimension - 1;

                        //使用横向的四点,求出X偏移量为amountX的临时插值点
                        value0 = GetValue(sourceVerticies, sourceDimension, iMinus1, indexJ, 0f);
                        value1 = GetValue(sourceVerticies, sourceDimension, i, indexJ, 1f);
                        value2 = GetValue(sourceVerticies, sourceDimension, iPlus1, indexJ, 2f);
                        value3 = GetValue(sourceVerticies, sourceDimension, iPlus2, indexJ, 3f);

                        tempValue[k] = MathHelper.CatmullRom(value0, value1, value2, value3, amountX);//插值
                    }

                    float vResult = MathHelper.CatmullRom(tempValue[0], tempValue[1], tempValue[2], tempValue[3], amountY);
                    destVerticies[y * destDimension + x] = new VertexPositionColorNormal(new Vector3((float)x - destDimension / 2, (float)y - destDimension / 2, vResult), Color.White, Vector3.Up);
                }
            }
            return destVerticies;
        }

        /// <summary>
        /// 灰度图像插值增加网络密度
        /// 这里先将源网格的高度值value转换为灰度bitmap,然后利用GDI的图片缩放功能放大图片到目标密度后,提取插值后的各像素点值作为新网格的各点高度值.
        /// </summary>
        /// <param name="sourceVerticies"></param>
        /// <param name="sourceDimension"></param>
        /// <param name="destDimension"></param>
        /// <returns></returns>
        public static VertexPositionColorNormal[] PopulateMeshGrayScaleBitmap(VertexPositionColorNormal[] sourceVerticies, int sourceDimension, int destDimension)
        {
            float min = float.MaxValue;
            float max = float.MinValue;
            for (int i = 0; i < sourceVerticies.Length; i++)
            {
                if (sourceVerticies[i].Position.Z < min) min = sourceVerticies[i].Position.Z;
                if (sourceVerticies[i].Position.Z > max) max = sourceVerticies[i].Position.Z;
            }
            //将源网格转换为位图
            Bitmap bitmapSource = new Bitmap(sourceDimension, sourceDimension, PixelFormat.Format24bppRgb);
            for (int x = 0; x < sourceDimension; x++)
            {
                for (int y = 0; y < sourceDimension; y++)
                {
                    float value = sourceVerticies[x * sourceDimension + y].Position.Z;
                    ColorD color = ValueToColor(min, max, value);
                    bitmapSource.SetPixel(x, y, color);
                }
            }

            Bitmap destBitmap = new Bitmap(destDimension, destDimension);
            Graphics g = Graphics.FromImage(destBitmap);

            // 插值算法的质量
            g.InterpolationMode = InterpolationMode.Bilinear;
            g.DrawImage(bitmapSource, new System.Drawing.Rectangle(0, 0, destDimension, destDimension), new System.Drawing.Rectangle(0, 0, bitmapSource.Width, bitmapSource.Height), GraphicsUnit.Pixel);
            g.Dispose();

            //输出测试图片
            //bitmapSource.Save("c:\\1.bmp");
            //destBitmap.Save("c:\\2.bmp");

            //将放大的图片还原为网格
            VertexPositionColorNormal[] destVerticies = new VertexPositionColorNormal[destDimension * destDimension];
            for (int x = 0; x < destDimension; x++)
            {
                for (int y = 0; y < destDimension; y++)
                {
                    float newValue = ColorToValue(destBitmap.GetPixel(x, y), min, max);
                    destVerticies[x * destDimension + y] = new VertexPositionColorNormal(new Vector3((float)x - destDimension / 2, (float)y - destDimension / 2, newValue), Microsoft.Xna.Framework.Graphics.Color.White, Vector3.Up);
                }
            }
            return destVerticies;
        }

        static ColorD ValueToColor(float min, float max, float value)
        {
            float amount = (value - min) / (max - min);
            Vector3 v0 = new Vector3(0, 0, 0);
            Vector3 v1 = new Vector3(255, 255, 255);
            Vector3 v = Vector3.Lerp(v0, v1, amount);
            ColorD color = ColorD.FromArgb((int)v.X, (int)v.Y, (int)v.Z);

            return color;
        }

        static float ColorToValue(ColorD color, float min, float max)
        {
            Vector3 v0 = new Vector3(0, 0, 0);
            Vector3 v1 = new Vector3(255, 255, 255);
            Vector3 v = new Vector3(color.R, color.G, color.B);
            Vector3 vectorTotal = (v1 - v0);
            return min + (max - min) * v.Length() / v1.Length();
        }

        static float GetValue(VertexPositionColorNormal[] sourceVerticies, int dimension, int i, int j, float x)
        {
            return sourceVerticies[j * dimension + i].Position.Z;
        }




        /// <summary>
        /// 生成点阵网格的三角形索引
        /// </summary>
        /// <param name="width"></param>
        /// <param name="height"></param>
        /// <returns></returns>
        public static ushort[] CreateIndex()
        {
            List<ushort> listIndex = new List<ushort>();
            for (ushort i = 0; i < GridDimension - 1; i++)
            {
                for (short j = 0; j < GridDimension - 1; j++)
                {
                    listIndex.Add((ushort)(i * GridDimension + j));
                    listIndex.Add((ushort)(i * GridDimension + j + 1));
                    listIndex.Add((ushort)((i + 1) * GridDimension + j));
                    listIndex.Add((ushort)((i + 1) * GridDimension + j));
                    listIndex.Add((ushort)(i * GridDimension + j + 1));
                    listIndex.Add((ushort)((i + 1) * GridDimension + j + 1));
                }
            }
            ushort[] arrIndex = listIndex.ToArray();
            return arrIndex;
        }

        ///// <summary>
        ///// 创建圆形的三角形索引,以矩形中心为圆心,以width/2为半径的圆
        ///// </summary>
        ///// <param name="width"></param>
        ///// <param name="height"></param>
        ///// <returns></returns>
        //public static int[] CreateCircleIndex(int width, int height)
        //{
        //    float radius = width / 2;
        //    int midX = width / 2;
        //    int midY = height / 2;
        //    List<int> listIndex = new List<int>();

        //    int[] arrStartEdegeIndex = new int[height];
        //    int[] arrEndEdegeIndex = new int[height];

        //    for (int i = 0; i < height - 1; i++)
        //    {
        //        bool inCircle = false;
        //        for (int j = 0; j < width - 1; j++)
        //        {
        //            int deltaX = j - midX;
        //            int deltaY = i - midY;
        //            if (Math.Sqrt(deltaX * deltaX + deltaY * deltaY) < radius) //在圆内
        //            {
        //                listIndex.Add(i * width + j);
        //                listIndex.Add(i * width + j + 1);
        //                listIndex.Add((i + 1) * width + j);
        //                listIndex.Add((i + 1) * width + j);
        //                listIndex.Add(i * width + j + 1);
        //                listIndex.Add((i + 1) * width + j + 1);
        //                if (inCircle == false) //刚进入圆
        //                {
        //                    inCircle = true;
        //                    arrStartEdegeIndex[i] = j;
        //                    //取第一个点
        //                }

        //            }
        //            else
        //            {
        //                if (inCircle == true)
        //                {
        //                    inCircle = false;
        //                    arrEndEdegeIndex[i] = j;
        //                }
        //            }
        //        }
        //    }


        //    for (int i = 0; i < height - 1; i++)
        //    {
        //        if (i < midY)
        //        {

        //            listIndex.Add(i * width + arrStartEdegeIndex[i]);
        //            listIndex.Add((i + 1) * width + arrStartEdegeIndex[i + 1]);
        //            listIndex.Add((i + 1) * width + arrStartEdegeIndex[i]);

        //            listIndex.Add(i * width + arrEndEdegeIndex[i] + 1);
        //            listIndex.Add((i + 1) * width + arrEndEdegeIndex[i] + 1);
        //            listIndex.Add((i + 1) * width + arrEndEdegeIndex[i + 1]);
        //        }
        //        else
        //        {
        //        }

        //    }


        //    int[] arrIndex = listIndex.ToArray();
        //    return arrIndex;
        //}



    }


}
