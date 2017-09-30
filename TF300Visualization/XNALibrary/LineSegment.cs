using System;
using System.Collections.Generic;
using System.Text;
using Microsoft.Xna.Framework;
using XNAHelper.Helpers;

namespace TF300.App.GUI.DatabaseUI.XNALibrary
{

    
    public class LineSegment
    {
        public Vector3[] Vetex = new Vector3[2];
         public LineSegment(Vector3 v0, Vector3 v1)
         {
             Vetex[0] = v0;
             Vetex[1] = v1;
         }

        /// <summary>
        /// 判断线段是否在一个平面的内侧
        /// </summary>
        /// <param name="plane"></param>
        /// <returns></returns>
        public bool OutsidePlane(Plane plane)
        {
            foreach (Vector3 vetex in Vetex)
            {
                if (PositionOutsidePlane(vetex, plane))
                    return true;
            }
            return false;

        }

        /// <summary>
        /// 判断点在平面的内侧还是外侧
        /// </summary>
        /// <param name="position"></param>
        /// <param name="plane"></param>
        /// <returns></returns>
        public static bool PositionOutsidePlane(Vector3 position, Plane plane)
        {
            Ray ray = new Ray(position, plane.Normal);
            float? d = ray.Intersects(plane);
            if (d.HasValue)//有交点则说明在外侧
            {
                if (d.Value > 0.05)
                    return true;
                else
                    return false;
            }
            else
            {
                return false;
            }

        }

        /// <summary>
        /// 用平面将线段分割为两段
        /// </summary>
        /// <param name="plane">平面</param>
        /// <returns>如果平面与线段相交则返回两段被平面分割的线段,否则返回原线段</returns>
        public List<LineSegment> SpliteByPlane(Plane plane)
        {
            List<LineSegment> listResult = new List<LineSegment>();
            Vector3 direction = Vetex[1] - Vetex[0];
            
            Ray ray = new Ray(Vetex[0], direction);
            float? d = ray.Intersects(plane);
            if (d.HasValue && d.Value<direction.Length())//射线有交点并且交点在线段以内
            {
                Vector3 intersectPoint = Vetex[0]+ Vector3.Normalize(direction) * d.Value;
                listResult.Add(new LineSegment(Vetex[0], intersectPoint));
                listResult.Add(new LineSegment(intersectPoint,Vetex[1]));
                return listResult;
            }
            listResult.Add(this);
            return listResult;
        }
    }
}
