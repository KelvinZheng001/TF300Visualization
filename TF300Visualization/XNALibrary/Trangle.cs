using System;
using System.Collections.Generic;
using System.Text;
using XNAHelper.CustomVertexs;
using Microsoft.Xna.Framework;
using XNAHelper.Helpers;

namespace TF300.App.GUI.DatabaseUI.XNALibrary
{
    public class Trangle
    {
        public Trangle()
        {
        }
        public Trangle(Vector3 v0, Vector3 v1, Vector3 v2)
        {
            //将三角形修正为同一时钟方向
            Vector3 firstVec = v0 - v1;
            Vector3 secondVec = v0 - v2;
            Vector3 normal = Vector3.Cross(secondVec, firstVec);
            if (normal.Z >= 0)
            {
                Vetex[0] = v0;
                Vetex[1] = v2;
                Vetex[2] = v1;
            }
            else
            {
                Vetex[0] = v0;
                Vetex[1] = v1;
                Vetex[2] = v2;
            }
        }
        public Vector3[] Vetex = new Vector3[3];



        /// <summary>
        /// 判断三角形与平面是否有交点
        /// </summary>
        /// <returns></returns>
        public bool IntersectsAltitude(float altitude)
        {
            bool AbovePoint = false;
            bool UnderPoint = false;
            for (int i = 0; i < 3; i++)
            {
                if (Vetex[i].Z < altitude)
                    UnderPoint = true;
                else
                    AbovePoint = true;

            }
            return UnderPoint && AbovePoint;
        }

        /// <summary>
        /// 返回三角形重心在XY平面上的投影到圆心的距离
        /// </summary>
        /// <returns></returns>
        public float DistanceXY()
        {
            Vector3 gravitipoint = Vector3.Zero;
            for (int i = 0; i < 3; i++)
            {
                gravitipoint += this.Vetex[i];
            }
            gravitipoint /= 3f;
            gravitipoint.Z = 0;
            return gravitipoint.Length();
        }

        ///// <summary>
        ///// 判断点在平面的内侧还是外侧
        ///// </summary>
        ///// <param name="position"></param>
        ///// <param name="plane"></param>
        ///// <returns></returns>
        //private bool PositionInsidePlane(Vector3 position,Plane plane)
        //{
        //    Ray ray = new Ray(position, plane.Normal);
        //    float? d = ray.Intersects(plane);
        //    if (d.HasValue )//有交点则说明在外侧
        //    {
        //        return false;
        //    }
        //    else
        //    {
        //        return true;
        //    }
            
        //}

        /// <summary>
        /// 判断三角形是否在平面内侧
        /// </summary>
        /// <param name="plane"></param>
        /// <returns></returns>
        public bool InsidePlane(Plane plane)
        {
            foreach (Vector3 vetex in Vetex)
            {
                if (PrimitiveHelper.PositionInsidePlane(vetex, plane))
                    return true;
            }
            return false;
        }


    //    /// <summary>
    //    /// 使用平面分割三角形
    //    /// </summary>
    //    /// <param name="plane"></param>
    //    public List<Trangle> DevideByPlane(Plane plane, ref List<Trangle> listTrangleInside, ref List<Trangle> listTrangleOutside, ref LineSegment IntersectLine)
    //    {
    //        List<Trangle> listResult = new List<Trangle>();
    //        listTrangleInside = new List<Trangle>();
    //        listTrangleOutside = new List<Trangle>();
    //        listResult.Add(this);

    //        Vector3[] V = new Vector3[3];
    //        for (int i = 0; i < 3; i++)
    //        {
    //            V[i] = this.Vetex[i];
    //        }

    //        List<int> listIndex = new List<int>();//相交边的序号

    //        List<Vector3> listStartV = new List<Vector3>();
    //        List<Vector3> listEndV = new List<Vector3>();

    //        //List<Vector3> listIntersectPoint = new List<Vector3>();
    //        Vector3[] arrIntersectPoint = new Vector3[3];

    //        listStartV.Add(V[0]);
    //        listStartV.Add(V[1]);
    //        listStartV.Add(V[2]);

    //        listEndV.Add(V[1]);
    //        listEndV.Add(V[2]);
    //        listEndV.Add(V[0]);
    //        int nCount = 0;
    //        for (int i = 0; i < 3; i++)
    //        {
    //            Vector3 V1 = listStartV[i];
    //            Vector3 V2 = listEndV[i];

    //            Vector3 direction = V2 - V1;
    //            float length = direction.Length();
    //            Vector3 directionNormalized = Vector3.Normalize(direction);


    //            Ray ray = new Ray(V1, directionNormalized);
    //            float? d = ray.Intersects(plane);
    //            if (d.HasValue && d > 0f)
    //            {
    //                if (d < length)//有交点
    //                {
    //                    Vector3 intersectPoint = V1 + directionNormalized * d.Value;
    //                    arrIntersectPoint[i] = intersectPoint;

    //                    nCount++;
    //                    listIndex.Add(i); //记住是哪条边相交
    //                }
    //            }
    //        }

    //        //List<Vector3> listButtonVector = new List<Vector3>();
    //        //Vector3 topVector;
    //        if (nCount != 0 && nCount != 2)
    //        {
    //            int i = nCount;
    //        }

    //        if (nCount == 2)//有相交
    //        {

    //            listResult.Clear();
    //            if (listIndex[0] == 0 && listIndex[1] == 1)
    //            {
    //                listResult.Add(new Trangle(Vetex[1], arrIntersectPoint[0], arrIntersectPoint[1]));
    //                listResult.Add(new Trangle(Vetex[0], arrIntersectPoint[0], arrIntersectPoint[1]));
    //                listResult.Add(new Trangle(Vetex[0], Vetex[2], arrIntersectPoint[1]));
    //                IntersectLine = new LineSegment(arrIntersectPoint[0], arrIntersectPoint[1]);
    //            }
    //            else if (listIndex[0] == 1 && listIndex[1] == 2)
    //            {
    //                listResult.Add(new Trangle(Vetex[2], arrIntersectPoint[1], arrIntersectPoint[2]));
    //                listResult.Add(new Trangle(Vetex[1], arrIntersectPoint[1], arrIntersectPoint[2]));
    //                listResult.Add(new Trangle(Vetex[0], Vetex[1], arrIntersectPoint[2]));
    //                IntersectLine = new LineSegment(arrIntersectPoint[1], arrIntersectPoint[2]);
    //            }
    //            else if (listIndex[0] == 0 && listIndex[1] == 2)
    //            {
    //                listResult.Add(new Trangle(Vetex[0], arrIntersectPoint[0], arrIntersectPoint[2]));
    //                listResult.Add(new Trangle(Vetex[1], arrIntersectPoint[2], arrIntersectPoint[0]));
    //                listResult.Add(new Trangle(Vetex[1], Vetex[2], arrIntersectPoint[2]));
    //                IntersectLine = new LineSegment(arrIntersectPoint[0], arrIntersectPoint[2]);
    //            }
    //            else
    //            {
    //                throw new Exception("trangle error");
    //            }

    //        }

    //        //计算处于平面内部的三角形
    //        foreach (Trangle trangle in listResult)
    //        {
    //            bool bInside = false;
    //            foreach (Vector3 vetex in trangle.Vetex)
    //            {
    //                if (PrimitiveHelper.PositionInsidePlane(vetex, plane))
    //                {
    //                    listTrangleInside.Add(trangle);
    //                    bInside = true;
    //                    break;
    //                }
    //            }
    //            if (!bInside)
    //            {
    //                listTrangleOutside.Add(trangle);
    //            }
    //        }
    //        return listResult;
    //    }

    //}

        /// <summary>
        /// 使用平面分割三角形
        /// </summary>
        /// <param name="plane"></param>
        public List<Trangle> SplitByPlane(Plane plane, ref LineSegment IntersectLine)
        {
            List<Trangle> listResult = new List<Trangle>();
            listResult.Add(this);

            Vector3[] V = new Vector3[3];
            for (int i = 0; i < 3; i++)
            {
                V[i] = this.Vetex[i];
            }

            List<int> listIndex = new List<int>();//相交边的序号

            List<Vector3> listStartV = new List<Vector3>();
            List<Vector3> listEndV = new List<Vector3>();

            //List<Vector3> listIntersectPoint = new List<Vector3>();
            Vector3[] arrIntersectPoint = new Vector3[3];

            listStartV.Add(V[0]);
            listStartV.Add(V[1]);
            listStartV.Add(V[2]);

            listEndV.Add(V[1]);
            listEndV.Add(V[2]);
            listEndV.Add(V[0]);
            int nCount = 0;
            for (int i = 0; i < 3; i++)
            {
                Vector3 V1 = listStartV[i];
                Vector3 V2 = listEndV[i];

                Vector3 direction = V2 - V1;
                float length = direction.Length();
                Vector3 directionNormalized = Vector3.Normalize(direction);


                Ray ray = new Ray(V1, directionNormalized);
                float? d = ray.Intersects(plane);
                if (d.HasValue && d > 0f)
                {
                    if (d < length)//有交点
                    {
                        Vector3 intersectPoint = V1 + directionNormalized * d.Value;
                        arrIntersectPoint[i] = intersectPoint;

                        nCount++;
                        listIndex.Add(i); //记住是哪条边相交
                    }
                }
            }

            //List<Vector3> listButtonVector = new List<Vector3>();
            //Vector3 topVector;
            if (nCount != 0 && nCount != 2)
            {
                int i = nCount;
            }

            if (nCount == 2)//有相交
            {

                listResult.Clear();
                if (listIndex[0] == 0 && listIndex[1] == 1)
                {
                    listResult.Add(new Trangle(Vetex[1], arrIntersectPoint[0], arrIntersectPoint[1]));
                    listResult.Add(new Trangle(Vetex[0], arrIntersectPoint[0], arrIntersectPoint[1]));
                    listResult.Add(new Trangle(Vetex[0], Vetex[2], arrIntersectPoint[1]));
                    IntersectLine = new LineSegment(arrIntersectPoint[0], arrIntersectPoint[1]);
                }
                else if (listIndex[0] == 1 && listIndex[1] == 2)
                {
                    listResult.Add(new Trangle(Vetex[2], arrIntersectPoint[1], arrIntersectPoint[2]));
                    listResult.Add(new Trangle(Vetex[1], arrIntersectPoint[1], arrIntersectPoint[2]));
                    listResult.Add(new Trangle(Vetex[0], Vetex[1], arrIntersectPoint[2]));
                    IntersectLine = new LineSegment(arrIntersectPoint[1], arrIntersectPoint[2]);
                }
                else if (listIndex[0] == 0 && listIndex[1] == 2)
                {
                    listResult.Add(new Trangle(Vetex[0], arrIntersectPoint[0], arrIntersectPoint[2]));
                    listResult.Add(new Trangle(Vetex[1], arrIntersectPoint[2], arrIntersectPoint[0]));
                    listResult.Add(new Trangle(Vetex[1], Vetex[2], arrIntersectPoint[2]));
                    IntersectLine = new LineSegment(arrIntersectPoint[0], arrIntersectPoint[2]);
                }
                else
                {
                    throw new Exception("trangle error");
                }

            }

            ////计算处于平面内部的三角形
            //foreach (Trangle trangle in listResult)
            //{
            //    bool bInside = false;
            //    foreach (Vector3 vetex in trangle.Vetex)
            //    {
            //        if (PrimitiveHelper.PositionInsidePlane(vetex, plane))
            //        {
            //            listTrangleInside.Add(trangle);
            //            bInside = true;
            //            break;
            //        }
            //    }
            //    if (!bInside)
            //    {
            //        listTrangleOutside.Add(trangle);
            //    }
            //}
            return listResult;
        }

    }
}
