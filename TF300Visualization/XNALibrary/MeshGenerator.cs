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
    //    /// ͨ��GridValue����
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
    //        //��������������
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
        /// ������������
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
            //���в�ֵ
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
        /// ת�����ݵ������Ϊ3Dͼ������
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
        ///// ������������,
        ///// </summary>
        ///// <param name="listPointValue">����������</param>
        ///// <returns>�����ʺ�XNA��ʾ�����ݸ�ʽ</returns>
        //public static CustomMesh GetMesh(List<PointValue> listPointValue)
        //{
        //    GridValue grid = GetGridData(listPointValue);
        //    CustomMesh mesh = new CustomMesh(grid);
        //    return mesh;
        //}


        /// <summary>
        /// �ÿ���𷨲�ֵ�������񶥵�
        /// </summary>
        /// <param name="listPointValue"></param>
        /// <param name="radius"></param>
        /// <returns></returns>
        public static VertexPositionColorNormal[] InterportationKrigingVertexPositionColor(List<PointValue> listPointValue, float radius)
        {


            GridValue grid = InterpolateToGridData(listPointValue, radius);//��ֵ��������
            VertexPositionColorNormal[] vertices = new VertexPositionColorNormal[grid.collumns * grid.rows];


            for (int i = 0; i < grid.collumns; i++)
            {
                for (int j = 0; j < grid.rows; j++)
                {
                    PointValue point = grid.ValueList[i * grid.collumns + j];
                    vertices[i * grid.collumns + j] = new VertexPositionColorNormal(new Vector3(point.X, point.Y, (float)point.Value), Color.Red, Vector3.Zero);
                }
            }

            //���������ܶ�
            //��ԭʼɢ�������Ƚ��нϵ��ܶ�InterpolateDimension�Ŀ�����ֵ�γ����������,��CatmullRom��ֵ�����������ܶ�,��ʹ���������ͼ��Ч��ƽ��
            VertexPositionColorNormal[] vertices1 = PopulateMeshCatmullRom(vertices, InterpolateDimension, GridDimension);
            return vertices1;


            //return vertices;
        }

        /// <summary>
        /// �������Ӽ��ܶ� ʹ��CatmullRom��ֵ
        /// </summary>
        /// <param name="sourceVerticies">Դ����</param>
        /// <param name="sourceDimension">Դ������ܶ�</param>
        /// <param name="destDimension">Ŀ���ܶ�</param>
        /// <returns></returns>
        public static VertexPositionColorNormal[] PopulateMeshCatmullRom(VertexPositionColorNormal[] sourceVerticies, int sourceDimension, int destDimension)
        {
            VertexPositionColorNormal[] destVerticies = new VertexPositionColorNormal[destDimension * destDimension];

            //x,y��Ŀ�������±�
            for (int y = 0; y < destDimension; y++)
            {
                for (int x = 0; x < destDimension; x++)
                {
                    //��Դ�����ȡֵ��
                    int i = x * sourceDimension / destDimension;
                    int j = y * sourceDimension / destDimension;
                    //���ֵ����x,y�����ϵ�ƫ����
                    float amountX = ((float)x / (float)destDimension - (float)i / (float)sourceDimension) / (1 / (float)sourceDimension);
                    float amountY = ((float)y / (float)destDimension - (float)j / (float)sourceDimension) / (1 / (float)sourceDimension);

                    float[] tempValue = new float[4];

                    float value0;
                    float value1;
                    float value2;
                    float value3;


                    for (int k = 0; k < 4; k++) //��ƫ����ΪamountX�������ĸ���ʱ��ֵ���ֵ
                    {
                        int iMinus1 = i - 1 < 0 ? i : i - 1;
                        int iPlus1 = i + 1 > sourceDimension - 1 ? i : i + 1;
                        int iPlus2 = iPlus1 + 1 > sourceDimension - 1 ? iPlus1 : iPlus1 + 1;

                        int indexJ = j - 1 + k;
                        if (indexJ < 0)
                            indexJ = 0;
                        if (indexJ >= sourceDimension)
                            indexJ = sourceDimension - 1;

                        //ʹ�ú�����ĵ�,���Xƫ����ΪamountX����ʱ��ֵ��
                        value0 = GetValue(sourceVerticies, sourceDimension, iMinus1, indexJ, 0f);
                        value1 = GetValue(sourceVerticies, sourceDimension, i, indexJ, 1f);
                        value2 = GetValue(sourceVerticies, sourceDimension, iPlus1, indexJ, 2f);
                        value3 = GetValue(sourceVerticies, sourceDimension, iPlus2, indexJ, 3f);

                        tempValue[k] = MathHelper.CatmullRom(value0, value1, value2, value3, amountX);//��ֵ
                    }

                    float vResult = MathHelper.CatmullRom(tempValue[0], tempValue[1], tempValue[2], tempValue[3], amountY);
                    destVerticies[y * destDimension + x] = new VertexPositionColorNormal(new Vector3((float)x - destDimension / 2, (float)y - destDimension / 2, vResult), Color.White, Vector3.Up);
                }
            }
            return destVerticies;
        }

        /// <summary>
        /// �Ҷ�ͼ���ֵ���������ܶ�
        /// �����Ƚ�Դ����ĸ߶�ֵvalueת��Ϊ�Ҷ�bitmap,Ȼ������GDI��ͼƬ���Ź��ܷŴ�ͼƬ��Ŀ���ܶȺ�,��ȡ��ֵ��ĸ����ص�ֵ��Ϊ������ĸ���߶�ֵ.
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
            //��Դ����ת��Ϊλͼ
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

            // ��ֵ�㷨������
            g.InterpolationMode = InterpolationMode.Bilinear;
            g.DrawImage(bitmapSource, new System.Drawing.Rectangle(0, 0, destDimension, destDimension), new System.Drawing.Rectangle(0, 0, bitmapSource.Width, bitmapSource.Height), GraphicsUnit.Pixel);
            g.Dispose();

            //�������ͼƬ
            //bitmapSource.Save("c:\\1.bmp");
            //destBitmap.Save("c:\\2.bmp");

            //���Ŵ��ͼƬ��ԭΪ����
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
        /// ���ɵ������������������
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
        ///// ����Բ�ε�����������,�Ծ�������ΪԲ��,��width/2Ϊ�뾶��Բ
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
        //            if (Math.Sqrt(deltaX * deltaX + deltaY * deltaY) < radius) //��Բ��
        //            {
        //                listIndex.Add(i * width + j);
        //                listIndex.Add(i * width + j + 1);
        //                listIndex.Add((i + 1) * width + j);
        //                listIndex.Add((i + 1) * width + j);
        //                listIndex.Add(i * width + j + 1);
        //                listIndex.Add((i + 1) * width + j + 1);
        //                if (inCircle == false) //�ս���Բ
        //                {
        //                    inCircle = true;
        //                    arrStartEdegeIndex[i] = j;
        //                    //ȡ��һ����
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
