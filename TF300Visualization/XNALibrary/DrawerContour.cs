using System;
using System.Collections.Generic;
using System.Text;
using Microsoft.Xna.Framework;
using XNAHelper.CustomVertexs;
using Microsoft.Xna.Framework.Graphics;
using XNAHelper.Helpers;
using System.Drawing;

namespace TF300.App.GUI.DatabaseUI.XNALibrary
{
    using Color = Microsoft.Xna.Framework.Graphics.Color;
    using ColorD = System.Drawing.Color;
    using TF300.App.GUI.DatabaseUI.XNALibrary.Helpers;
    using Microsoft.Xna.Framework.Content;
    public class DrawerContourMap : Drawer
    {


        bool _Flat;
        /// <summary>
        /// ָ���Ƿ����Ϊƽ��
        /// </summary>
        public bool Flat
        {
            get { return _Flat; }
            set
            {
                _Flat = value;
                _OptionApplyed = false;
            }
        }

        bool _DrawContourLine;
        /// <summary>
        /// ָ���Ƿ���Ƶȸ���
        /// </summary>
        public bool DrawContourLine
        {
            get { return _DrawContourLine; }
            set
            {
                _DrawContourLine = value;
                _OptionApplyed = false;
            }
        }

        bool _DrawContourZone;
        /// <summary>
        /// ָ���Ƿ���Ƶȸ�����
        /// </summary>
        public bool DrawContourZone
        {
            get { return _DrawContourZone; }
            set
            {
                _DrawContourZone = value;
                _OptionApplyed = false;
            }
        }


        VertexPositionColorNormal[] mVertices;
        VertexPositionColorNormal[] mEdgeVertices;
        ushort[] mIndices;
        List<Color> mListColor = new List<Color>();

        ColorPalette mColorPalette = new ColorPalette();
        int _LayerCount = 10;
        /// <summary>
        /// ָ���ȸ��滮�ֲ���
        /// </summary>
        public int LayerCount
        {
            get { return _LayerCount; }
            set { _LayerCount = value; }
        }
        float mMin;
        float mMax;
        float mDelta;

        VertexPositionColorNormal[] mContourVertices;
        ushort[] mContourIndices;
        ushort[] mEdgeIndices;



        //��ͼ����

        VertexBuffer vertexBufferContour;
        IndexBuffer indexBufferContour;
        VertexBuffer vertexBufferBody;
        IndexBuffer indexBufferBody;
        VertexBuffer vertexBufferEdge;
        IndexBuffer indexBufferEdge;







        public DrawerContourMap(GraphicsDevice graphicsDevice, ContentManager Content, VertexPositionColorNormal[] verticies, ushort[] indices, VertexPositionColorNormal[] edgeVerticies,int nLayerCount)
            : base(graphicsDevice, Content)
        {

            //mListColor.Add(Color.Red);
            //mListColor.Add(Color.Tomato);
            //mListColor.Add(Color.Wheat);
            //mListColor.Add(Color.SteelBlue);
            //mListColor.Add(Color.SkyBlue);
            //mListColor.Add(Color.LightYellow);
            //mListColor.Add(Color.Yellow);
            //mListColor.Add(Color.YellowGreen);
            //mListColor.Add(Color.MediumSeaGreen);
            //mListColor.Add(Color.SpringGreen);
            //mListColor.Add(Color.SteelBlue);
            //mListColor.Add(Color.Purple);
            //mListColor.Add(Color.Plum);
            //mListColor.Add(Color.PeachPuff);
            //mListColor.Add(Color.PaleTurquoise);
            //mListColor.Add(Color.Orange);
            //mListColor.Add(Color.Navy);
            //mListColor.Add(Color.MintCream);
            //mListColor.Add(Color.MediumBlue);
            //mListColor.Add(Color.LimeGreen);
            //mListColor.Add(Color.LightSkyBlue);

            //mColorPalette
            _LayerCount = nLayerCount;
            for (int i = 0; i < _LayerCount; i++)
            {
                mListColor.Add(mColorPalette.GetColor((float)i/(float)(_LayerCount-1)));
            }
                

            //_LayerCount = mListColor.Count - 1;

            //���ƶ�����������
            mVertices = new VertexPositionColorNormal[verticies.Length];
            verticies.CopyTo(mVertices, 0);
            mIndices = indices;
            mEdgeVertices = new VertexPositionColorNormal[edgeVerticies.Length];
            edgeVerticies.CopyTo(mEdgeVertices, 0);

            mEdgeIndices = new ushort[edgeVerticies.Length];
            for (ushort i = 0; i < mEdgeIndices.Length; i++)
            {
                mEdgeIndices[i] = i;
            }


            float minValue = float.MaxValue;
            float maxValue = float.MinValue;
            foreach (VertexPositionColorNormal Vertex in verticies)
            {
                if (minValue > Vertex.Position.Z) minValue = Vertex.Position.Z;
                if (maxValue < Vertex.Position.Z) maxValue = Vertex.Position.Z;
            }
            mMin = minValue;
            mMax = maxValue;
            mDelta = maxValue - minValue;

            //����ͼ��
            InitializeLegend();

            //����ȸ���
            CalculateContourLine();
            Calculate();

            //mDrawVertices = new VertexPositionColorNormal[mVertices.Length];
            //mVertices.CopyTo(mDrawVertices, 0);

            //���ŵ�һ����Χ
            PrimitiveHelper.ScaleZToRange(mVertices, 10f, minValue, maxValue);

            //���㷨��
            PrimitiveHelper.GenerateNormalsForTriangleList(mVertices, mIndices);

            PrimitiveHelper.ScaleZToRange(mContourVertices, 10f, minValue, maxValue);
            PrimitiveHelper.ScaleZToRange(mEdgeVertices, 10f, minValue, maxValue);

            //InitializeEffect(graphicsDevice);
        }

        List<string> mListAltitudeString = new List<string>();
        Vector2 mLegendPosition = new Vector2(20, 120);

        public Vector2 LegendPosition
        {
            get { return mLegendPosition; }
            set { mLegendPosition = value; }
        }
        ///// <summary>
        ///// ����ͼ��
        ///// </summary>
        //private void CalculateLegend()
        //{
        //    int layerCount = mListColor.Count - 1;

        //}



        Dictionary<float, VertexPositionColorNormal[]> mDictContourVertexByAltitude = new Dictionary<float, VertexPositionColorNormal[]>();
        Dictionary<float, int[]> mDictContourIndexByAltitude = new Dictionary<float, int[]>();
        /// <summary>
        /// ����ȸ���
        /// </summary>
        private void CalculateContourLine()
        {

            //float testMin = float.MaxValue;
            //float testMax = float.MinValue;
            //foreach (VertexPositionColorNormal Vertex in mVertices)
            //{
            //    if (testMin > Vertex.Position.Z) testMin = Vertex.Position.Z;
            //    if (testMax < Vertex.Position.Z) testMax = Vertex.Position.Z;
            //}

            float fAltitude = 0f; //�߶�

            for (int layer = 0; layer <= _LayerCount; layer++)//ѭ��ÿ���ȸ߲�
            {

                fAltitude = (float)layer * mDelta / (float)_LayerCount + mMin;
                mListAltitudeString.Add(fAltitude.ToString("f5"));


                //����ȸ���

                int nTrangleCount = mIndices.Length / 3;//����������
                List<Vector3> listAboveVertex = new List<Vector3>();
                List<Vector3> listUnderVertex = new List<Vector3>();
                Plane plane = new Plane(Vector3.Forward, fAltitude);
                List<Vector3> listContourPoint = new List<Vector3>();
                //����ÿһ��������
                for (int i = 0; i < nTrangleCount; i++)
                {
                    listAboveVertex.Clear();
                    listUnderVertex.Clear();
                    //�������Ϊ����ƽ�������ƽ���������
                    for (int j = 0; j < 3; j++)
                    {
                        VertexPositionColorNormal Vertex = mVertices[mIndices[i * 3 + j]];
                        if (Vertex.Position.Z >= fAltitude)
                        {
                            listAboveVertex.Add(Vertex.Position);
                        }
                        else
                        {
                            listUnderVertex.Add(Vertex.Position);
                        }
                    }

                    //Ray testRay = new Ray(new Vector3(0, 0, 16), Vector3.Forward);
                    //Plane testPlane = new Plane(Vector3.Backward, -15);
                    //float? d1 = testRay.Intersects(testPlane);
                    //Vector3 point = Vector3.Zero + Vector3.Forward * d1.Value;

                    //����Ϸ��е�,�·�Ҳ�е�,��˵��ƽ�����������ཻ.
                    if (listAboveVertex.Count > 0 && listUnderVertex.Count > 0)
                    {
                        List<Vector3> listIntersectPoint = new List<Vector3>();//�����б�
                        //���н������
                        foreach (Vector3 v1 in listAboveVertex)
                        {
                            foreach (Vector3 v2 in listUnderVertex)
                            {
                                Vector3 direction = Vector3.Normalize(v2 - v1);
                                Ray ray = new Ray(v1, direction);
                                float? d = ray.Intersects(plane);
                                if (d.HasValue)
                                {
                                    Vector3 Point = v1 + direction * d.Value;
                                    listIntersectPoint.Add(Point);
                                }
                            }
                        }
                        if (listIntersectPoint.Count == 2)
                        {
                            listContourPoint.Add(listIntersectPoint[0]);
                            listContourPoint.Add(listIntersectPoint[1]);
                        }
                    }
                }
                VertexPositionColorNormal[] arrVertexContour = new VertexPositionColorNormal[listContourPoint.Count];
                int[] arrIndexContour = new int[listContourPoint.Count];
                for (int i = 0; i < listContourPoint.Count; i++)
                {
                    Color color = GetColor(fAltitude);
                    arrVertexContour[i] = new VertexPositionColorNormal(listContourPoint[i], color, Vector3.Up);
                    arrIndexContour[i] = i;//���ɸ߶���������
                }
                if (arrVertexContour.Length > 0)
                {
                    mDictContourVertexByAltitude.Add(fAltitude, arrVertexContour);
                    mDictContourIndexByAltitude.Add(fAltitude, arrIndexContour);
                }
            }//�߶�ѭ��

            int nVerticesCount = 0;
            foreach (VertexPositionColorNormal[] vertices in mDictContourVertexByAltitude.Values)
            {
                nVerticesCount += vertices.Length;
            }//�����ܳ�
            mContourVertices = new VertexPositionColorNormal[nVerticesCount];
            mContourIndices = new ushort[nVerticesCount];
            ushort index = 0;
            foreach (VertexPositionColorNormal[] vertices in mDictContourVertexByAltitude.Values)
            {
                for (int i = 0; i < vertices.Length; i++)
                {
                    mContourVertices[index] = vertices[i];
                    mContourIndices[index] = index;
                    index++;
                }
            }



        }

        Texture2D mLegendTexture;
        int mLegendWidth = 30;
        int mLegendHeight = 20;

        private void InitializeLegend()
        {
            //����ͼ������ͼ����
            System.Drawing.Bitmap bmp = new System.Drawing.Bitmap(
                                               mLegendWidth, mLegendHeight * _LayerCount,
                                               System.Drawing.Imaging.PixelFormat.Format24bppRgb
                                             );
            Graphics graphics = Graphics.FromImage(bmp);
            SolidBrush brush = new SolidBrush(ColorD.Black);
            for (int i = 0; i < _LayerCount; i++)
            {
                ColorD color = ColorD.FromArgb(mListColor[i].R, mListColor[i].G, mListColor[i].B);
                brush.Color = color;
                graphics.FillRectangle(brush, new System.Drawing.Rectangle(0, mLegendHeight * (_LayerCount-i-1), mLegendWidth, mLegendHeight));
            }
            brush.Dispose();
            graphics.Dispose();
            //λ��ת��������
            mLegendTexture = TextureHelper.Texture2DFromBitmap(bmp, this.mGraphicsDevice);

        }

        /// <summary>
        /// ����ȸ���
        /// </summary>
        private void Calculate()
        {

            List<VertexPositionColorNormal> listVertices = new List<VertexPositionColorNormal>(mVertices);
            List<VertexPositionColorNormal> listGerneratedVertices = new List<VertexPositionColorNormal>();
            List<int> listGerneratedIndices = new List<int>();
            List<Trangle> ListNonIntersectTrangles = new List<Trangle>();
            List<Trangle> ListIntersectTrangles = new List<Trangle>();

            List<Trangle> listTrangle = new List<Trangle>();
            int nCount = mIndices.Length / 3;//����������
            //����ÿһ��������,ת��Ϊ������trangle��������.
            for (int i = 0; i < nCount; i++)
            {
                Trangle trangle = new Trangle(mVertices[mIndices[i * 3]].Position, mVertices[mIndices[i * 3 + 1]].Position, mVertices[mIndices[i * 3 + 2]].Position);
                listTrangle.Add(trangle);
            }



            //����ÿ��������,��ÿ���������ж������и߶Ȳ��Ƿ��ཻ,��������ཻ,��������ཻ������. �ⲽ����Ҳ��,��Ϊ�˵�һ�ξ�ȥ��������������
            foreach (Trangle trangle in listTrangle)
            {

                for (int layer = 0; layer <= _LayerCount; layer++)//ѭ��ÿ���ȸ߲�
                {
                    float fAltitude = (float)layer * mDelta / (float)_LayerCount + mMin;
                    List<Vector3> listAboveVertex = new List<Vector3>();
                    List<Vector3> listUnderVertex = new List<Vector3>();
                    //�������Ϊ����ƽ�������ƽ���������
                    for (int j = 0; j < 3; j++)
                    {
                        Vector3 vertex = trangle.Vetex[j];
                        if (vertex.Z >= fAltitude)
                        {
                            listAboveVertex.Add(vertex);
                        }
                        else
                        {
                            listUnderVertex.Add(vertex);
                        }
                    }
                    if (listAboveVertex.Count > 0 && listUnderVertex.Count > 0)
                    {
                        ListIntersectTrangles.Add(trangle);
                        break;
                    }
                }//foreach layer
            }//foreach trangle 
            //ʣ�ಿ�ּ������ཻ�������б���
            foreach (Trangle trangle in listTrangle)
            {
                if (!ListIntersectTrangles.Contains(trangle))
                {
                    ListNonIntersectTrangles.Add(trangle);
                }
            }

            listTrangle.Clear();//������������б����ڴ�ż���������.


            //�����ཻ������

            for (int layer = 0; layer <= _LayerCount; layer++)//ѭ��ÿ���ȸ߲�
            {
                List<Trangle> listTrangleTemp = new List<Trangle>();//��ʱ�����δ����
                float fAltitude = (float)layer * mDelta / (float)_LayerCount + mMin;
                float fNextAltitude = (float)(layer + 1) * mDelta / (float)_LayerCount + mMin;
                Plane plane = new Plane(Vector3.Forward, fAltitude);
                foreach (Trangle trangle in ListIntersectTrangles)
                {
                    List<Vector3> listAboveVertex = new List<Vector3>();
                    List<Vector3> listUnderVertex = new List<Vector3>();

                    //�������Ϊ����ƽ�������ƽ���������
                    for (int j = 0; j < 3; j++)
                    {
                        Vector3 vertex = trangle.Vetex[j];
                        if (vertex.Z >= fAltitude)
                        {
                            listAboveVertex.Add(vertex);
                        }
                        else
                        {
                            listUnderVertex.Add(vertex);
                        }
                    }
                    if (listAboveVertex.Count > 0 && listUnderVertex.Count > 0)
                    {
                        //Vector3[] IntersectPoint = new Vector3[2];
                        Vector3[] AbovePoint = new Vector3[2];
                        Vector3[] UnderPoint = new Vector3[2];

                        int index = 0;
                        List<Vector3> listIntersectPoint = new List<Vector3>();//�����б�
                        //���н������
                        foreach (Vector3 AboveVertex in listAboveVertex)
                        {
                            foreach (Vector3 UnderVertex in listUnderVertex)
                            {
                                Vector3 direction = Vector3.Normalize(UnderVertex - AboveVertex);
                                Ray ray = new Ray(AboveVertex, direction);
                                float? d = ray.Intersects(plane);
                                if (d.HasValue)
                                {
                                    Vector3 Point = AboveVertex + direction * d.Value;
                                    listIntersectPoint.Add(Point);
                                    AbovePoint[index] = AboveVertex;
                                    UnderPoint[index] = UnderVertex;//���������һ�ε�������Ϊ�˱�����IntersectPoint��������˳��һ��
                                    index++;
                                }

                            }
                        }
                        if (listIntersectPoint.Count == 2)
                        {
                            //������������·��㹹�ɵ���������������б�,���Ϸ����������㹹�ɵ��������������ʱ�洢��

                            if (listUnderVertex.Count == 1)
                            {
                                //�·�
                                Trangle newTrangle = new Trangle(
                                    listIntersectPoint[0],
                                    listIntersectPoint[1],
                                    UnderPoint[0]
                                );
                                listTrangle.Add(newTrangle);
                                //�Ϸ�������ʱ����
                                newTrangle = new Trangle(listIntersectPoint[0], AbovePoint[0], AbovePoint[1]
                                );
                                if (newTrangle.IntersectsAltitude(fNextAltitude))
                                {
                                    listTrangleTemp.Add(newTrangle);
                                }
                                else
                                {
                                    listTrangle.Add(newTrangle);
                                }
                                newTrangle = new Trangle(listIntersectPoint[1], listIntersectPoint[0], AbovePoint[1]
                                );
                                if (newTrangle.IntersectsAltitude(fNextAltitude))
                                {
                                    listTrangleTemp.Add(newTrangle);
                                }
                                else
                                {
                                    listTrangle.Add(newTrangle);
                                }

                            }
                            else
                            {
                                //�·�
                                Trangle newTrangle = new Trangle(listIntersectPoint[0], listIntersectPoint[1], UnderPoint[0]
                                );
                                listTrangle.Add(newTrangle);
                                newTrangle = new Trangle(listIntersectPoint[1], UnderPoint[0], UnderPoint[1]
                                );
                                listTrangle.Add(newTrangle);
                                //�Ϸ�
                                newTrangle = new Trangle(listIntersectPoint[1], listIntersectPoint[0], AbovePoint[0]
                                );
                                if (newTrangle.IntersectsAltitude(fNextAltitude))
                                {
                                    listTrangleTemp.Add(newTrangle);
                                }
                                else
                                {
                                    listTrangle.Add(newTrangle);
                                }

                            }
                        }
                    }
                    else //�뵱ǰ�߶�û�н������������ʱ�洢�б�,������һ�ּ���
                    {
                        listTrangleTemp.Add(trangle);
                    }
                }//foreach trangle
                //ʹ����ʱ�洢�����������б�����һ������
                ListIntersectTrangles = listTrangleTemp;
            }//foreach layer


            //�����ཻ�������μ��뵽���б�
            foreach (Trangle trangle in ListNonIntersectTrangles)
            {
                listTrangle.Add(trangle);
            }


            //�������δ��б�ת��Ϊ��������
            mVertices = new VertexPositionColorNormal[listTrangle.Count * 3];
            for (int i = 0; i < listTrangle.Count; i++)
            {
                Trangle trangle = listTrangle[i];
                for (int j = 0; j < 3; j++)
                {
                    //���������ε���ɫ: ����ȡ�����ε�һ����,�ж�������һ��
                    float value = (trangle.Vetex[0].Z + trangle.Vetex[1].Z + trangle.Vetex[2].Z) / 3;
                    Color color = GetColor(value);
                    mVertices[i * 3 + j] = new VertexPositionColorNormal(trangle.Vetex[j], color, Vector3.Up);
                }
            }
            mIndices = new ushort[listTrangle.Count * 3];
            for (ushort i = 0; i < mIndices.Length; i++)
            {
                if (i == ushort.MaxValue) break;
                mIndices[i] = i;
            }


        }

        /// <summary>
        /// ����ĳ���߶ȵ���ɫ
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        private Color GetColor(float value)
        {

            if (value < mMin)
            {
                value = mMin;
            }
            int layerCount = _LayerCount;
            int layer = 0;
            if (mDelta == 0)
            {
                layerCount = 0;
            }
            else
            {
                layer = (int)((value - mMin) / ((mDelta) / layerCount));
                if (layer > mListColor.Count - 1)
                    layer = mListColor.Count - 1;
            }
            return mListColor[layer];
        }


        VertexPositionColorNormal[] mDrawVertices;
        VertexPositionColorNormal[] mDrawContourVertices;
        VertexPositionColorNormal[] mDrawEdgeVertices;


        protected override void ClearBuffer()
        {

            if (indexBufferContour != null)
            {
                indexBufferContour.Dispose();
            }

            if (vertexBufferContour != null)
            {
                vertexBufferContour.Dispose();
            }

            if (vertexBufferBody != null)
            {
                vertexBufferBody.Dispose();
            }
            if (indexBufferBody != null)
            {
                indexBufferBody.Dispose();
            }

            if (vertexBufferEdge != null)
                vertexBufferEdge.Dispose();
            if (indexBufferEdge != null)
                indexBufferEdge.Dispose();

            base.ClearBuffer();
        }



        /// <summary>
        /// ������Ӧ�õ�ͼ��
        /// </summary>
        protected override void ApplyOptions()
        {

            ClearBuffer();

            mDrawVertices = new VertexPositionColorNormal[mVertices.Length];
            mVertices.CopyTo(mDrawVertices, 0);
            mDrawContourVertices = new VertexPositionColorNormal[mContourVertices.Length];
            mContourVertices.CopyTo(mDrawContourVertices, 0);
            mDrawEdgeVertices = new VertexPositionColorNormal[mEdgeVertices.Length];
            mEdgeVertices.CopyTo(mDrawEdgeVertices, 0);

            if (_Flat)//ƽ����ʾ
            {
                for (int i = 0; i < mDrawVertices.Length; i++)
                {
                    mDrawVertices[i].Position.Z = 0;
                    mDrawVertices[i].Normal = Vector3.Up;
                }
                for (int i = 0; i < mDrawContourVertices.Length; i++)
                {
                    mDrawContourVertices[i].Position.Z = 0;
                    mDrawContourVertices[i].Normal = Vector3.Up;
                }
                for (int i = 0; i < mDrawEdgeVertices.Length; i++)
                {
                    mDrawEdgeVertices[i].Position.Z = 0;
                }
            }



            // Create a vertex declaration, describing the format of our vertex data.
            vertexDeclaration = new VertexDeclaration(mGraphicsDevice, VertexPositionColorNormal.VertexElements);

            if (_DrawContourLine)
            {
                if (this.mDrawContourVertices.Length > 0)
                {
                    //�ȸ��߶��㻺��
                    vertexBufferContour = new VertexBuffer(mGraphicsDevice, typeof(VertexPositionColorNormal), this.mDrawContourVertices.Length, BufferUsage.None);
                    vertexBufferContour.SetData(mDrawContourVertices);
                    indexBufferContour = new IndexBuffer(mGraphicsDevice, typeof(ushort), this.mContourIndices.Length, BufferUsage.None);
                    indexBufferContour.SetData(mContourIndices);
                    //��Ե��
                    vertexBufferEdge = new VertexBuffer(mGraphicsDevice, typeof(VertexPositionColorNormal), mDrawEdgeVertices.Length, BufferUsage.None);
                    vertexBufferEdge.SetData(mDrawEdgeVertices);
                    indexBufferEdge = new IndexBuffer(mGraphicsDevice, typeof(ushort), this.mEdgeIndices.Length, BufferUsage.None);
                    indexBufferEdge.SetData(mEdgeIndices);
                }

            }

            if (_DrawContourZone)
            {
                //ģ���嶥�㻺��
                vertexBufferBody = new VertexBuffer(mGraphicsDevice, typeof(VertexPositionColorNormal), this.mDrawVertices.Length, BufferUsage.None);
                vertexBufferBody.SetData(mDrawVertices);
                indexBufferBody = new IndexBuffer(mGraphicsDevice, typeof(ushort), this.mIndices.Length, BufferUsage.None);
                indexBufferBody.SetData(mIndices);
            }

            _OptionApplyed = true;

        }


        /// <summary>
        /// ����ͼ��
        /// </summary>
        protected override void Draw(BasicEffect effect)
        {


            GraphicsDevice graphicsDevice = effect.GraphicsDevice;
            RenderState renderState = graphicsDevice.RenderState;
            renderState.FillMode = FillMode.Solid;
            renderState.CullMode = CullMode.None;

            //���ö�������������
            graphicsDevice.VertexDeclaration = vertexDeclaration;
            //graphicsDevice.Vertices[0].SetSource(vertexBufferContour, 0, VertexPositionColorNormal.SizeInBytes);
            //graphicsDevice.Indices = indexBufferContour;




            //Vector3 point4 = graphicsDevice.Viewport.Unproject(new Vector3(location.X, location.Y + 5, 0f), effect.Projection, effect.View, effect.World);





            effect.Begin();

            foreach (EffectPass effectPass in effect.CurrentTechnique.Passes)
            {
                effectPass.Begin();

                //�ȸ���
                if (_DrawContourLine && mDrawContourVertices.Length > 0)
                {
                    //�ȸ���
                    graphicsDevice.Vertices[0].SetSource(vertexBufferContour, 0, VertexPositionColorNormal.SizeInBytes);
                    graphicsDevice.Indices = indexBufferContour;
                    graphicsDevice.DrawIndexedPrimitives(PrimitiveType.LineList, 0, 0, mContourVertices.Length, 0, mContourIndices.Length / 2);
                    //��Ե��
                    graphicsDevice.Vertices[0].SetSource(vertexBufferEdge, 0, VertexPositionColorNormal.SizeInBytes);
                    graphicsDevice.Indices = indexBufferEdge;
                    graphicsDevice.DrawIndexedPrimitives(PrimitiveType.LineList, 0, 0, this.mEdgeVertices.Length, 0, this.mEdgeIndices.Length / 2);


                }
                if (_DrawContourZone)//ģ���� �ȸ�����
                {
                    graphicsDevice.Vertices[0].SetSource(vertexBufferBody, 0, VertexPositionColorNormal.SizeInBytes);
                    graphicsDevice.Indices = indexBufferBody;
                    graphicsDevice.DrawIndexedPrimitives(PrimitiveType.TriangleList, 0, 0, this.mDrawVertices.Length, 0, mIndices.Length / 3);
                }


                effectPass.End();
            }

            //����ͼ��
            mSpriteBatch.Begin();

            mSpriteBatch.Draw(mLegendTexture, mLegendPosition, Color.White);
            for (int i = 0; i < _LayerCount; i++)
            {
                mSpriteBatch.DrawString(mFont, mListAltitudeString[this._LayerCount-i], new Vector2(mLegendPosition.X + mLegendWidth, mLegendPosition.Y + mLegendHeight * i),Color.SteelBlue, 0, new Vector2(0, 0), 0.7f, SpriteEffects.None, 0);
            }

            mSpriteBatch.End();


            effect.End();
        }


    }
}
