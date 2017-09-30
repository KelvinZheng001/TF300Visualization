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
        /// �ж��߶��Ƿ���һ��ƽ����ڲ�
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
        /// �жϵ���ƽ����ڲ໹�����
        /// </summary>
        /// <param name="position"></param>
        /// <param name="plane"></param>
        /// <returns></returns>
        public static bool PositionOutsidePlane(Vector3 position, Plane plane)
        {
            Ray ray = new Ray(position, plane.Normal);
            float? d = ray.Intersects(plane);
            if (d.HasValue)//�н�����˵�������
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
        /// ��ƽ�潫�߶ηָ�Ϊ����
        /// </summary>
        /// <param name="plane">ƽ��</param>
        /// <returns>���ƽ�����߶��ཻ�򷵻����α�ƽ��ָ���߶�,���򷵻�ԭ�߶�</returns>
        public List<LineSegment> SpliteByPlane(Plane plane)
        {
            List<LineSegment> listResult = new List<LineSegment>();
            Vector3 direction = Vetex[1] - Vetex[0];
            
            Ray ray = new Ray(Vetex[0], direction);
            float? d = ray.Intersects(plane);
            if (d.HasValue && d.Value<direction.Length())//�����н��㲢�ҽ������߶�����
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
