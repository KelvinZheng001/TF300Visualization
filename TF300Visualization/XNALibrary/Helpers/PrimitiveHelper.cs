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
        /// ���㶥��ķ���
        /// </summary>
        /// <param name="vertices"></param>
        /// <param name="indices"></param>
        /// <returns>���㷨������</returns>
        public static Vector3[] GenerateNormalsForTriangleList(VertexPositionColor[] vertices, ushort[] indices)
        {
            Vector3[] vertexNormals = new Vector3[vertices.Length];
            //��ʼ����������
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
        /// ���㶥��ķ���
        /// </summary>
        /// <param name="vertices"></param>
        /// <param name="indices"></param>
        /// <returns></returns>
        public static void GenerateNormalsForTriangleList(VertexPositionColorNormal[] vertices, ushort[] indices)
        {
            Dictionary<Vector3, Vector3> dictNormalByPosition = new Dictionary<Vector3, Vector3>();
            //��ʼ����������
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
        //    //��ʼ����������
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
        /// Ϊģ���Ͻ���ɫ
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


            float deltaHeight = maxValue - minValue;//��С�����Ĳ�ֵ

            for (int i = 0; i < vertices.Length; i++)
            {
                VertexPositionColorNormal vertex = vertices[i];
                double value = vertex.Position.Z;
                //������ɫ ��ɫ����Ϊ������ʼɫ����ɫ. ����Ϊ��ֵ��߶�ֵ������Ȩֵ.
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
        /// ��Z�������ŵ�ָ����Χ
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
                ZScale = ZRange / (maxValue - minValue); //����Ŵ���
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
        /// �Զ������Բ���μ���
        /// </summary>
        public static void ClipByCylinder(ref VertexPositionColorNormal[] vertices, ref ushort[] indices, int gridDimension, ref VertexPositionColorNormal[] edgeVertices)
        {
            float radius = ((float)gridDimension - 2.2f) / 2f;

            List<Trangle> listTrangle = new List<Trangle>();
            List<Trangle> listTrangleAdd = new List<Trangle>();//����Ѽ��úõ�������
            List<Trangle> listTrangleRemove = new List<Trangle>();

            List<LineSegment> listEdgeLine = new List<LineSegment>();

            List<Trangle> List1 = new List<Trangle>();
            List<Trangle> List2 = new List<Trangle>();

            int nCount = indices.Length / 3;//����������
            //����ÿһ��������,ת��Ϊ������trangle��������.
            for (int i = 0; i < nCount; i++)
            {
                Trangle trangle = new Trangle(vertices[indices[i * 3]].Position, vertices[indices[i * 3 + 1]].Position, vertices[indices[i * 3 + 2]].Position);
                listTrangle.Add(trangle);
            }

            int tessellation = 33;

            //�ָ������ཻ��������




            for (int i = 0; i < tessellation; i++)//ÿ����
            {
                Vector3 normal = -GetCircleVector(i, tessellation);
                //normal.X = -normal.X;
                //normal.Y = -normal.Y;
                //normal.Z = -normal.Z;
                Plane plane = new Plane(normal, radius);


                //�����Ե��
                for (int edgelineIndex = listEdgeLine.Count - 1; edgelineIndex >= 0; edgelineIndex--)
                {
                    LineSegment segment = listEdgeLine[edgelineIndex];
                    List<LineSegment> listSegment = segment.SpliteByPlane(plane);
                    if (listSegment.Count == 2)
                    {
                        listEdgeLine.RemoveAt(edgelineIndex);//�Ƴ�ԭ�߶�
                        foreach (LineSegment splitSegment in listSegment)
                        {
                            if (!splitSegment.OutsidePlane(plane))//��ƽ������
                            {
                                listEdgeLine.Add(splitSegment);//������ƽ���ڵ��߶�
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
                    if (listDividedTrangle.Count == 1)//û���ཻ
                    {
                        if (!trangle.InsidePlane(plane)) //��ƽ�����,��ɾ��
                        {
                            listTrangle.RemoveAt(trangleIndex);
                        }
                    }
                    else if (listDividedTrangle.Count == 3)//�ཻ���ѱ����
                    {
                        foreach (Trangle devideTrangle in listDividedTrangle)
                        {
                            if (devideTrangle.InsidePlane(plane))
                            {
                                listTrangle.Add(devideTrangle);
                            }
                        }
                        listEdgeLine.Add(edgeLine);//�ཻ�������ε��ཻ����Ϊ��Ե��,�����б�
                        listTrangle.RemoveAt(trangleIndex);//�Ƴ�ԭ������
                    }
                }
            }

            for (int i = 0; i < tessellation; i++)//ÿ����
            {
                Vector3 normal = -GetCircleVector(i, tessellation);
                Plane plane = new Plane(normal, radius);

                //�����Ե��
                for (int edgelineIndex = listEdgeLine.Count - 1; edgelineIndex >= 0; edgelineIndex--)
                {
                    LineSegment segment = listEdgeLine[edgelineIndex];
                    if (segment.OutsidePlane(plane))//��ƽ������,ɾ��
                    {
                        listEdgeLine.RemoveAt(edgelineIndex);
                    }
                }
            }



            //�������δ��б�ת��Ϊ��������
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


            //����Ե���б�ת��Ϊ��������
            edgeVertices = new VertexPositionColorNormal[listEdgeLine.Count * 2];
            for (int i = 0; i < listEdgeLine.Count; i++)
            {
                LineSegment segment = listEdgeLine[i];
                edgeVertices[i * 2] = new VertexPositionColorNormal(segment.Vetex[0], Color.Black, Vector3.Up);
                edgeVertices[i * 2 + 1] = new VertexPositionColorNormal(segment.Vetex[1], Color.Black, Vector3.Up);
            }

        }


        /// <summary>
        /// ȡԲ�������ĸ�������
        /// </summary>
        /// <param name="i">Ҫȡ�Ķκ�</param>
        /// <param name="tessellation">�ָ����</param>
        /// <returns></returns>
        static Vector3 GetCircleVector(int i, int tessellation)
        {
            float angle = i * MathHelper.TwoPi / tessellation + 0.1f;
            float dx = (float)Math.Cos(angle);
            float dy = (float)Math.Sin(angle);

            return new Vector3(dx, dy, 0);
        }



        /// <summary>
        /// �жϵ���ƽ����ڲ໹�����
        /// </summary>
        /// <param name="position"></param>
        /// <param name="plane"></param>
        /// <returns></returns>
        public static bool PositionInsidePlane(Vector3 position, Plane plane)
        {
            Ray ray = new Ray(position, plane.Normal);
            float? d = ray.Intersects(plane);
            if (d.HasValue)//�н�����˵�������
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
