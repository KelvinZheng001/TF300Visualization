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
        /// 指定是否绘制为平面
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
        /// 指定是否绘制等高线
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
        /// 指定是否绘制等高区域
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
        /// 指定等高面划分层数
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



        //绘图所需

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

            //复制顶点数组数据
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

            //计算图例
            InitializeLegend();

            //计算等高线
            CalculateContourLine();
            Calculate();

            //mDrawVertices = new VertexPositionColorNormal[mVertices.Length];
            //mVertices.CopyTo(mDrawVertices, 0);

            //缩放到一定范围
            PrimitiveHelper.ScaleZToRange(mVertices, 10f, minValue, maxValue);

            //计算法线
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
        ///// 计算图例
        ///// </summary>
        //private void CalculateLegend()
        //{
        //    int layerCount = mListColor.Count - 1;

        //}



        Dictionary<float, VertexPositionColorNormal[]> mDictContourVertexByAltitude = new Dictionary<float, VertexPositionColorNormal[]>();
        Dictionary<float, int[]> mDictContourIndexByAltitude = new Dictionary<float, int[]>();
        /// <summary>
        /// 计算等高线
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

            float fAltitude = 0f; //高度

            for (int layer = 0; layer <= _LayerCount; layer++)//循环每个等高层
            {

                fAltitude = (float)layer * mDelta / (float)_LayerCount + mMin;
                mListAltitudeString.Add(fAltitude.ToString("f5"));


                //计算等高线

                int nTrangleCount = mIndices.Length / 3;//三角形数量
                List<Vector3> listAboveVertex = new List<Vector3>();
                List<Vector3> listUnderVertex = new List<Vector3>();
                Plane plane = new Plane(Vector3.Forward, fAltitude);
                List<Vector3> listContourPoint = new List<Vector3>();
                //遍历每一个三角型
                for (int i = 0; i < nTrangleCount; i++)
                {
                    listAboveVertex.Clear();
                    listUnderVertex.Clear();
                    //将顶点分为高于平面与低于平面的两部分
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

                    //如果上方有点,下方也有点,则说明平面与三角形相交.
                    if (listAboveVertex.Count > 0 && listUnderVertex.Count > 0)
                    {
                        List<Vector3> listIntersectPoint = new List<Vector3>();//交点列表
                        //进行交点计算
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
                    arrIndexContour[i] = i;//生成高度连线索引
                }
                if (arrVertexContour.Length > 0)
                {
                    mDictContourVertexByAltitude.Add(fAltitude, arrVertexContour);
                    mDictContourIndexByAltitude.Add(fAltitude, arrIndexContour);
                }
            }//高度循环

            int nVerticesCount = 0;
            foreach (VertexPositionColorNormal[] vertices in mDictContourVertexByAltitude.Values)
            {
                nVerticesCount += vertices.Length;
            }//计算总长
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
            //生成图例的贴图纹理
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
            //位置转换到纹理
            mLegendTexture = TextureHelper.Texture2DFromBitmap(bmp, this.mGraphicsDevice);

        }

        /// <summary>
        /// 计算等高线
        /// </summary>
        private void Calculate()
        {

            List<VertexPositionColorNormal> listVertices = new List<VertexPositionColorNormal>(mVertices);
            List<VertexPositionColorNormal> listGerneratedVertices = new List<VertexPositionColorNormal>();
            List<int> listGerneratedIndices = new List<int>();
            List<Trangle> ListNonIntersectTrangles = new List<Trangle>();
            List<Trangle> ListIntersectTrangles = new List<Trangle>();

            List<Trangle> listTrangle = new List<Trangle>();
            int nCount = mIndices.Length / 3;//三角形数量
            //遍历每一个三角型,转换为独立的trangle对象数组.
            for (int i = 0; i < nCount; i++)
            {
                Trangle trangle = new Trangle(mVertices[mIndices[i * 3]].Position, mVertices[mIndices[i * 3 + 1]].Position, mVertices[mIndices[i * 3 + 2]].Position);
                listTrangle.Add(trangle);
            }



            //遍历每个三角形,对每个三角型判断与所有高度层是否相交,如果都不相交,则加入无相交三角型. 这步不做也行,是为了第一次就去掉多数的三角型
            foreach (Trangle trangle in listTrangle)
            {

                for (int layer = 0; layer <= _LayerCount; layer++)//循环每个等高层
                {
                    float fAltitude = (float)layer * mDelta / (float)_LayerCount + mMin;
                    List<Vector3> listAboveVertex = new List<Vector3>();
                    List<Vector3> listUnderVertex = new List<Vector3>();
                    //将顶点分为高于平面与低于平面的两部分
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
            //剩余部分加入无相交三角形列表中
            foreach (Trangle trangle in listTrangle)
            {
                if (!ListIntersectTrangles.Contains(trangle))
                {
                    ListNonIntersectTrangles.Add(trangle);
                }
            }

            listTrangle.Clear();//清空主三角形列表用于存放计算后的数据.


            //计算相交三角形

            for (int layer = 0; layer <= _LayerCount; layer++)//循环每个等高层
            {
                List<Trangle> listTrangleTemp = new List<Trangle>();//临时三角形存放区
                float fAltitude = (float)layer * mDelta / (float)_LayerCount + mMin;
                float fNextAltitude = (float)(layer + 1) * mDelta / (float)_LayerCount + mMin;
                Plane plane = new Plane(Vector3.Forward, fAltitude);
                foreach (Trangle trangle in ListIntersectTrangles)
                {
                    List<Vector3> listAboveVertex = new List<Vector3>();
                    List<Vector3> listUnderVertex = new List<Vector3>();

                    //将顶点分为高于平面与低于平面的两部分
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
                        List<Vector3> listIntersectPoint = new List<Vector3>();//交点列表
                        //进行交点计算
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
                                    UnderPoint[index] = UnderVertex;//这里再添加一次到数组是为了保持与IntersectPoint的两个点顺序一致
                                    index++;
                                }

                            }
                        }
                        if (listIntersectPoint.Count == 2)
                        {
                            //由两个交与点下方点构成的三角形输出到大列表,由上方点与两交点构成的三角形输出到临时存储区

                            if (listUnderVertex.Count == 1)
                            {
                                //下方
                                Trangle newTrangle = new Trangle(
                                    listIntersectPoint[0],
                                    listIntersectPoint[1],
                                    UnderPoint[0]
                                );
                                listTrangle.Add(newTrangle);
                                //上方加入临时三角
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
                                //下方
                                Trangle newTrangle = new Trangle(listIntersectPoint[0], listIntersectPoint[1], UnderPoint[0]
                                );
                                listTrangle.Add(newTrangle);
                                newTrangle = new Trangle(listIntersectPoint[1], UnderPoint[0], UnderPoint[1]
                                );
                                listTrangle.Add(newTrangle);
                                //上方
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
                    else //与当前高度没有交点则输出到临时存储列表,参与下一轮计算
                    {
                        listTrangleTemp.Add(trangle);
                    }
                }//foreach trangle
                //使用临时存储区的三角型列表作下一层运算
                ListIntersectTrangles = listTrangleTemp;
            }//foreach layer


            //将无相交的三角形加入到大列表
            foreach (Trangle trangle in ListNonIntersectTrangles)
            {
                listTrangle.Add(trangle);
            }


            //将三角形大列表转换为顶点数组
            mVertices = new VertexPositionColorNormal[listTrangle.Count * 3];
            for (int i = 0; i < listTrangle.Count; i++)
            {
                Trangle trangle = listTrangle[i];
                for (int j = 0; j < 3; j++)
                {
                    //计算三角形的颜色: 任意取三角形的一个点,判断属于哪一层
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
        /// 返回某个高度的颜色
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
        /// 将属性应用到图像
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

            if (_Flat)//平面显示
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
                    //等高线顶点缓冲
                    vertexBufferContour = new VertexBuffer(mGraphicsDevice, typeof(VertexPositionColorNormal), this.mDrawContourVertices.Length, BufferUsage.None);
                    vertexBufferContour.SetData(mDrawContourVertices);
                    indexBufferContour = new IndexBuffer(mGraphicsDevice, typeof(ushort), this.mContourIndices.Length, BufferUsage.None);
                    indexBufferContour.SetData(mContourIndices);
                    //边缘线
                    vertexBufferEdge = new VertexBuffer(mGraphicsDevice, typeof(VertexPositionColorNormal), mDrawEdgeVertices.Length, BufferUsage.None);
                    vertexBufferEdge.SetData(mDrawEdgeVertices);
                    indexBufferEdge = new IndexBuffer(mGraphicsDevice, typeof(ushort), this.mEdgeIndices.Length, BufferUsage.None);
                    indexBufferEdge.SetData(mEdgeIndices);
                }

            }

            if (_DrawContourZone)
            {
                //模形体顶点缓冲
                vertexBufferBody = new VertexBuffer(mGraphicsDevice, typeof(VertexPositionColorNormal), this.mDrawVertices.Length, BufferUsage.None);
                vertexBufferBody.SetData(mDrawVertices);
                indexBufferBody = new IndexBuffer(mGraphicsDevice, typeof(ushort), this.mIndices.Length, BufferUsage.None);
                indexBufferBody.SetData(mIndices);
            }

            _OptionApplyed = true;

        }


        /// <summary>
        /// 绘制图型
        /// </summary>
        protected override void Draw(BasicEffect effect)
        {


            GraphicsDevice graphicsDevice = effect.GraphicsDevice;
            RenderState renderState = graphicsDevice.RenderState;
            renderState.FillMode = FillMode.Solid;
            renderState.CullMode = CullMode.None;

            //设置顶点与索引缓存
            graphicsDevice.VertexDeclaration = vertexDeclaration;
            //graphicsDevice.Vertices[0].SetSource(vertexBufferContour, 0, VertexPositionColorNormal.SizeInBytes);
            //graphicsDevice.Indices = indexBufferContour;




            //Vector3 point4 = graphicsDevice.Viewport.Unproject(new Vector3(location.X, location.Y + 5, 0f), effect.Projection, effect.View, effect.World);





            effect.Begin();

            foreach (EffectPass effectPass in effect.CurrentTechnique.Passes)
            {
                effectPass.Begin();

                //等高线
                if (_DrawContourLine && mDrawContourVertices.Length > 0)
                {
                    //等高线
                    graphicsDevice.Vertices[0].SetSource(vertexBufferContour, 0, VertexPositionColorNormal.SizeInBytes);
                    graphicsDevice.Indices = indexBufferContour;
                    graphicsDevice.DrawIndexedPrimitives(PrimitiveType.LineList, 0, 0, mContourVertices.Length, 0, mContourIndices.Length / 2);
                    //边缘线
                    graphicsDevice.Vertices[0].SetSource(vertexBufferEdge, 0, VertexPositionColorNormal.SizeInBytes);
                    graphicsDevice.Indices = indexBufferEdge;
                    graphicsDevice.DrawIndexedPrimitives(PrimitiveType.LineList, 0, 0, this.mEdgeVertices.Length, 0, this.mEdgeIndices.Length / 2);


                }
                if (_DrawContourZone)//模形体 等高区域
                {
                    graphicsDevice.Vertices[0].SetSource(vertexBufferBody, 0, VertexPositionColorNormal.SizeInBytes);
                    graphicsDevice.Indices = indexBufferBody;
                    graphicsDevice.DrawIndexedPrimitives(PrimitiveType.TriangleList, 0, 0, this.mDrawVertices.Length, 0, mIndices.Length / 3);
                }


                effectPass.End();
            }

            //绘制图例
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
