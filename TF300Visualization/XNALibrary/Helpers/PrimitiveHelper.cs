using System;
using System.Collections.Generic;
using System.Text;
using Microsoft.Xna.Framework.Graphics;
using Microsoft.Xna.Framework;
using XNAHelper.CustomVertexs;
using TF300.App.GUI.DatabaseUI.XNALibrary;

namespace XNAHelper.Helpers
{
    public class PrimitiveHelper
    {
        /// <summary>
        /// 计算顶点的法线
        /// </summary>
        /// <param name="vertices"></param>
        /// <param name="indices"></param>
        /// <returns>顶点法线数组</returns>
        public static Vector3[] GenerateNormalsForTriangleList(VertexPositionColor[] vertices, ushort[] indices)
        {
            Vector3[] vertexNormals = new Vector3[vertices.Length];
            //初始化顶点数组
            for (int i = 0; i < vertices.Length; i++)
            {
                vertexNormals[i] = Vector3.Zero;
            }
            //vertices[i].Normal = new Vector3(0, 0, 0);

            for (int i = 0; i < indices.Length / 3; i++)
            {
                Vector3 firstVec = vertices[indices[i * 3 + 1]].Position - vertices[indices[i * 3]].Position;
                Vector3 secondVec = vertices[indices[i * 3 + 2]].Position - vertices[indices[i * 3]].Position;
                Vector3 normal = Vector3.Cross(secondVec, firstVec);
                normal.Normalize();

                vertexNormals[indices[i * 3]] += normal;
                vertexNormals[indices[i * 3 + 1]] += normal;
                vertexNormals[indices[i * 3 + 2]] += normal;
            }
            return vertexNormals;
        }//GenerateNormalsForTriangleList




        /// <summary>
        /// 计算顶点的法线
        /// </summary>
        /// <param name="vertices"></param>
        /// <param name="indices"></param>
        /// <returns></returns>
        public static void GenerateNormalsForTriangleList(VertexPositionColorNormal[] vertices, ushort[] indices)
        {
            Dictionary<Vector3, Vector3> dictNormalByPosition = new Dictionary<Vector3, Vector3>();
            //初始化顶点数组
            for (int i = 0; i < vertices.Length; i++)
            {
                vertices[i].Normal = Vector3.Zero;
            }
            //vertices[i].Normal = new Vector3(0, 0, 0);

            foreach (VertexPositionColorNormal vertext in vertices)
            {
                dictNormalByPosition[vertext.Position] = Vector3.Zero;
            }

            for (int i = 0; i < indices.Length / 3; i++)
            {
                Vector3 firstVec = vertices[indices[i * 3 + 1]].Position - vertices[indices[i * 3]].Position;
                Vector3 secondVec = vertices[indices[i * 3 + 2]].Position - vertices[indices[i * 3]].Position;
                Vector3 normal = Vector3.Cross(secondVec, firstVec);
                normal.Normalize();

                //vertices[indices[i * 3]].Normal += normal;
                //vertices[indices[i * 3 + 1]].Normal += normal;
                //vertices[indices[i * 3 + 2]].Normal += normal;

                dictNormalByPosition[vertices[indices[i * 3]].Position] = normal;
                dictNormalByPosition[vertices[indices[i * 3 + 1]].Position] = normal;
                dictNormalByPosition[vertices[indices[i * 3 + 2]].Position] = normal;

            }
            for (int i = 0; i < vertices.Length; i++)
            {
                vertices[i].Normal = dictNormalByPosition[vertices[i].Position];
            }

        }//GenerateNormalsForTriangleList

        //public static void GenerateNormalsForTriangleList(VertexPositionColorNormal[] vertices, int[] indices)
        //{
        //    //Vector3[] vertexNormals = new Vector3[vertices.Length];
        //    //初始化顶点数组
        //    for (int i = 0; i < vertices.Length; i++)
        //    {
        //        vertices[i].Normal = Vector3.Zero;
        //    }
        //    //vertices[i].Normal = new Vector3(0, 0, 0);

        //    for (int i = 0; i < indices.Length / 3; i++)
        //    {
        //        Vector3 firstVec = vertices[indices[i * 3 + 1]].Position - vertices[indices[i * 3]].Position;
        //        Vector3 secondVec = vertices[indices[i * 3 + 2]].Position - vertices[indices[i * 3]].Position;
        //        Vector3 normal = Vector3.Cross(secondVec, firstVec);
        //        normal.Normalize();

        //        vertices[indices[i * 3]].Normal += normal;
        //        vertices[indices[i * 3 + 1]].Normal += normal;
        //        vertices[indices[i * 3 + 2]].Normal += normal;
        //    }
        //}//GenerateNormalsForTriangleList
        /// <summary>
        /// 为模形上渐变色
        /// </summary>
        /// <param name="vertices"></param>
        /// <param name="indices"></param>
        public static void RenderGradientColor(VertexPositionColorNormal[] vertices, Color colorMiddle, Color colorTop, Color colorBotton)
        {


            float minValue = float.MaxValue;
            float maxValue = float.MinValue;
            foreach (VertexPositionColorNormal Vertex in vertices)
            {
                if (minValue > Vertex.Position.Z) minValue = Vertex.Position.Z;
                if (maxValue < Vertex.Position.Z) maxValue = Vertex.Position.Z;
            }


            float deltaHeight = maxValue - minValue;//最小与最大的差值

            for (int i = 0; i < vertices.Length; i++)
            {
                VertexPositionColorNormal vertex = vertices[i];
                double value = vertex.Position.Z;
                //计算颜色 颜色方向为等于起始色结束色. 距离为插值点高度值的线性权值.
                float distnace = (((float)value) - minValue) / deltaHeight;
                //Vector3 vectorColor = Vector3.CatmullRom(colorV1.ToVector3(), colorStart.ToVector3(), colorEnd.ToVector3(), colorV4.ToVector3(),-0.3f);
                Vector3 vectorColor = Vector3.Zero;
                if (distnace > 0.5)
                {
                    vectorColor = Vector3.Lerp(colorMiddle.ToVector3(), colorTop.ToVector3(), (distnace - 0.5f) * 2f);
                }
                else
                {
                    vectorColor = Vector3.Lerp(colorBotton.ToVector3(), colorMiddle.ToVector3(), distnace * 2f);
                }
                //Vector3 vectorColor = colorStart.ToVector3() + colorDirection * distnace;
                vertex.Color = new Color(vectorColor);//
                vertices[i] = vertex;
            }

        }

        /// <summary>
        /// 将Z坐标缩放到指定范围
        /// </summary>
        /// <param name="verticies"></param>
        /// <param name="ZRange"></param>
        /// <param name="maxValue"></param>
        /// <param name="minValue"></param>
        public static void ScaleZToRange(VertexPositionColorNormal[] verticies, float ZRange, float minValue, float maxValue)
        {

            //float minValue = float.MaxValue;
            //float maxValue = float.MinValue;
            //foreach (VertexPositionColorNormal Vertex in verticies)
            //{
            //    if (minValue > Vertex.Position.Z) minValue = Vertex.Position.Z;
            //    if (maxValue < Vertex.Position.Z) maxValue = Vertex.Position.Z;
            //}


            float ZScale;
            if (maxValue != minValue)
            {
                ZScale = ZRange / (maxValue - minValue); //计算放大倍数
            }
            else
            {
                ZScale = 0f;
            }

            for (int i = 0; i < verticies.Length; i++)
            {
                verticies[i].Position.Z = (verticies[i].Position.Z - minValue) * ZScale;
            }

        }



        /// <summary>
        /// 对顶点进行圆柱形剪裁
        /// </summary>
        public static void ClipByCylinder(ref VertexPositionColorNormal[] vertices, ref ushort[] indices, int gridDimension, ref VertexPositionColorNormal[] edgeVertices)
        {
            float radius = ((float)gridDimension - 2.2f) / 2f;

            List<Trangle> listTrangle = new List<Trangle>();
            List<Trangle> listTrangleAdd = new List<Trangle>();//存放已剪裁好的三角形
            List<Trangle> listTrangleRemove = new List<Trangle>();

            List<LineSegment> listEdgeLine = new List<LineSegment>();

            List<Trangle> List1 = new List<Trangle>();
            List<Trangle> List2 = new List<Trangle>();

            int nCount = indices.Length / 3;//三角形数量
            //遍历每一个三角型,转换为独立的trangle对象数组.
            for (int i = 0; i < nCount; i++)
            {
                Trangle trangle = new Trangle(vertices[indices[i * 3]].Position, vertices[indices[i * 3 + 1]].Position, vertices[indices[i * 3 + 2]].Position);
                listTrangle.Add(trangle);
            }

            int tessellation = 33;

            //分割所有相交的三角形




            for (int i = 0; i < tessellation; i++)//每个面
            {
                Vector3 normal = -GetCircleVector(i, tessellation);
                //normal.X = -normal.X;
                //normal.Y = -normal.Y;
                //normal.Z = -normal.Z;
                Plane plane = new Plane(normal, radius);


                //处理边缘线
                for (int edgelineIndex = listEdgeLine.Count - 1; edgelineIndex >= 0; edgelineIndex--)
                {
                    LineSegment segment = listEdgeLine[edgelineIndex];
                    List<LineSegment> listSegment = segment.SpliteByPlane(plane);
                    if (listSegment.Count == 2)
                    {
                        listEdgeLine.RemoveAt(edgelineIndex);//移除原线段
                        foreach (LineSegment splitSegment in listSegment)
                        {
                            if (!splitSegment.OutsidePlane(plane))//在平面以内
                            {
                                listEdgeLine.Add(splitSegment);//加入在平内内的线段
                            }
                        }
                    }
                }


                for (int trangleIndex = listTrangle.Count - 1; trangleIndex >= 0; trangleIndex--)
                {
                    //List<Trangle> listInsideTrangle = null;
                    //List<Trangle> listOutsideTrangle = null;
                    LineSegment edgeLine = null;
                    Trangle trangle = listTrangle[trangleIndex];
                    List<Trangle> listDividedTrangle = trangle.SplitByPlane(plane, ref edgeLine);
                    if (listDividedTrangle.Count == 1)//没有相交
                    {
                        if (!trangle.InsidePlane(plane)) //在平面外侧,则删除
                        {
                            listTrangle.RemoveAt(trangleIndex);
                        }
                    }
                    else if (listDividedTrangle.Count == 3)//相交并已被拆分
                    {
                        foreach (Trangle devideTrangle in listDividedTrangle)
                        {
                            if (devideTrangle.InsidePlane(plane))
                            {
                                listTrangle.Add(devideTrangle);
                            }
                        }
                        listEdgeLine.Add(edgeLine);//相交的三角形的相交线作为边缘线,加入列表
                        listTrangle.RemoveAt(trangleIndex);//移除原三角形
                    }
                }
            }

            for (int i = 0; i < tessellation; i++)//每个面
            {
                Vector3 normal = -GetCircleVector(i, tessellation);
                Plane plane = new Plane(normal, radius);

                //处理边缘线
                for (int edgelineIndex = listEdgeLine.Count - 1; edgelineIndex >= 0; edgelineIndex--)
                {
                    LineSegment segment = listEdgeLine[edgelineIndex];
                    if (segment.OutsidePlane(plane))//在平面以外,删除
                    {
                        listEdgeLine.RemoveAt(edgelineIndex);
                    }
                }
            }



            //将三角形大列表转换为顶点数组
            vertices = new VertexPositionColorNormal[listTrangle.Count * 3];
            for (int i = 0; i < listTrangle.Count; i++)
            {
                Trangle trangle = listTrangle[i];
                Color color = listTrangleAdd.Contains(trangle) ? Color.Red : Color.Green;
                for (int j = 0; j < 3; j++)
                {
                    vertices[i * 3 + j] = new VertexPositionColorNormal(trangle.Vetex[j], color, Vector3.Up);
                }
            }
            indices = new ushort[listTrangle.Count * 3];
            for (ushort i = 0; i < indices.Length; i++)
            {
                indices[i] = i;
            }


            //将边缘线列表转换为顶点数组
            edgeVertices = new VertexPositionColorNormal[listEdgeLine.Count * 2];
            for (int i = 0; i < listEdgeLine.Count; i++)
            {
                LineSegment segment = listEdgeLine[i];
                edgeVertices[i * 2] = new VertexPositionColorNormal(segment.Vetex[0], Color.Black, Vector3.Up);
                edgeVertices[i * 2 + 1] = new VertexPositionColorNormal(segment.Vetex[1], Color.Black, Vector3.Up);
            }

        }


        /// <summary>
        /// 取圆周向量的辅助方法
        /// </summary>
        /// <param name="i">要取的段号</param>
        /// <param name="tessellation">分割段数</param>
        /// <returns></returns>
        static Vector3 GetCircleVector(int i, int tessellation)
        {
            float angle = i * MathHelper.TwoPi / tessellation + 0.1f;
            float dx = (float)Math.Cos(angle);
            float dy = (float)Math.Sin(angle);

            return new Vector3(dx, dy, 0);
        }



        /// <summary>
        /// 判断点在平面的内侧还是外侧
        /// </summary>
        /// <param name="position"></param>
        /// <param name="plane"></param>
        /// <returns></returns>
        public static bool PositionInsidePlane(Vector3 position, Plane plane)
        {
            Ray ray = new Ray(position, plane.Normal);
            float? d = ray.Intersects(plane);
            if (d.HasValue)//有交点则说明在外侧
            {
                return false;
            }
            else
            {
                return true;
            }

        }

    }
}
